require 'rbbt/sources/organism'
module Immunomics

  input :mutated_isoforms, :array, "Mutated Isoforms"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :mutated_isoform_sequences => :tsv do |mis,organism|
    mi_sequence = Organism.protein_sequence(organism).tsv :persist => true
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Wildtype", "Pos", "WT", "Alt"], :type => :list, :organism => organism
    dumper.init
    TSV.traverse mis, :type => :array, :into => dumper do |mi|
      protein, _sep, change = mi.partition(":")
      wt, pos, alt = change.partition /\d+/
      wt_sequence = mi_sequence[protein]
      raise "No sequence for protein: #{ protein }" if wt_sequence.nil?
      raise "Wildtype aa at position #{pos} is #{wt_sequence[pos.to_i - 1]} not #{wt}" if wt != wt_sequence[pos.to_i - 1] 
      [mi, [wt_sequence, pos, wt, alt]]
    end
  end

  dep :mutated_isoform_sequences
  input :flank_size, :integer, "Flank size", 20
  task :mutation_flanking_sequence => :tsv do |flank_size|
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Lost", "Position of Change", "N Flank", "Change", "C Flank"], :type => :list
    dumper.init
    TSV.traverse step(:mutated_isoform_sequences), :into => dumper do |mi,values|
      mi = mi.first if Array === mi

      wt_sequence, pos, wt, alt = values
      pos = pos.to_i

      mut_sequence = wt_sequence.dup
      next if alt == "*"

      aa_size = 1
      lost = ""
      if m = alt.match(/(Indel|FrameShift)\((.*)\)/i)
        type = m[1]
        aas = m[2]

        next if aas.nil?
        next if aas == ""
        next if aas == "?"

        aas = "-" + aas unless %(+ -).include? aas[0]

        aas_a = aas.split("")
        while aas_a.first == '-'
          lost << mut_sequence[pos - 1]
          mut_sequence[pos - 1] = ""
          aas_a.shift
        end

        if aas_a.first == '+'
          pos += 1
          aas_a.shift
        end

        if type == "Indel"
          aa_size = aas_a.length
          mut_sequence[pos - 1] = (aas_a * "") + mut_sequence[pos - 1]  if aas_a
        else
          aa_size = 0
          lost = ""
          mut_sequence[pos - 1..-1] = aas_a * ""
        end
      else
        type = 'SNV'
        lost << mut_sequence[pos - 1]
        mut_sequence[pos - 1] = alt
      end

      n_start = pos-flank_size-1
      n_start = 0 if n_start < 0

      if pos > 1
        n_flank = mut_sequence[n_start..pos-2]
      else
        n_flank = ""
      end

      if type === "Indel" || aa_size > 0 # Indel
        pos_end = pos + aa_size - 1
        pos_str = mut_sequence[pos-1..pos_end-1]
        c_flank = mut_sequence[pos_end..pos+flank_size+aa_size-2]
      else # FrameShift
        pos_str = mut_sequence[pos-1..-1]
        pos_str << "-" unless pos_str[-1] == "*"
        c_flank = wt_sequence[(pos-1)..pos+flank_size-2]
      end

      n_flank = "-" * (flank_size - n_flank.length) + n_flank
      c_flank = c_flank + "-" * (flank_size - c_flank.length)
      [mi, [lost, pos, n_flank, pos_str, c_flank]]
    end
  end

  dep :mutation_flanking_sequence, :compute => :produce
  input :sizes, :array, "Peptide size", [9]
  task :epitopes => :tsv do |sizes|
    organism = self.recursive_inputs[:organism]
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Wildtype", "Mutated",  "Offset", "Lost", "Change", "N-term", "C-term"], :type => :double, :namespace => organism
    dumper.init
    TSV.traverse step(:mutation_flanking_sequence), :into => dumper do |mi,info|
      mi = mi.first if Array === mi
      lost, change_position, nterm, change, cterm = info

      frameshift = mi.include? "FrameShift"

      offset = nterm.length
      mut_peptide = nterm + change + cterm
      wt_peptide = nterm + lost + cterm

      epitopes = []
      offset = nterm.length
      diff = change.length - lost.length

      sizes.each do |size|
        size = size.to_i
        middle = change.length + (size - change.length) / 2 
        steps = size + change.length - 1
        steps.times do |step|
          index = offset + step - (size - 1)
          mut_epitope = mut_peptide[index..index+size - 1]
          next if mut_epitope.include? '-'
          next if mut_epitope.include? '*'

          left_wt_epitope = wt_peptide[index..index+size - 1]
          left_wt_epitope = "" if left_wt_epitope.nil? || left_wt_epitope.nil?
          right_wt_epitope = wt_peptide[index-diff..index+size-diff-1]
          right_wt_epitope = "" if frameshift || right_wt_epitope.nil?

          right_wt_epitope = "" if right_wt_epitope.include? "*"
          left_wt_epitope = "" if left_wt_epitope.include? "*"
          right_wt_epitope = "" if right_wt_epitope.length != size
          left_wt_epitope = "" if left_wt_epitope.length != size

          left_wt_epitope_common = 0
          size.times do |i|
            left_wt_epitope_common += 1 if left_wt_epitope[i] == mut_epitope[i]
          end

          right_wt_epitope_common = 0
          size.times do |i|
            right_wt_epitope_common += 1 if right_wt_epitope[i] == mut_epitope[i]
          end

          wt_epitope = if right_wt_epitope_common > left_wt_epitope_common
                         right_wt_epitope
                       else
                         left_wt_epitope
                       end

          wt_epitope = "" if [right_wt_epitope_common, left_wt_epitope_common].max <= [3, size / 2].max

          next if mut_epitope == wt_epitope

          epitopes << [wt_epitope, mut_epitope, step]
        end
      end

      [mi, Misc.zip_fields(epitopes) + [lost, change, nterm, cterm]]
    end
  end

end

