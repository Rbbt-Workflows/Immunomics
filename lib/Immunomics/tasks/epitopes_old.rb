module Immunomics
  dep :mutated_isoform_sequences
  input :sizes, :array, "Peptide size", [9]
  task :epitopes_old => :tsv do |sizes|
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Wildtype", "Mutated", "Offset", "Flank"], :type => :double
    dumper.init
    TSV.traverse step(:mutated_isoform_sequences), :into => dumper do |mi,values|
      mi = mi.first if Array === mi
      epitopes = []
      wt_sequence, pos, wt, alt = values
      mut_sequence = wt_sequence.dup
      next if alt == "*"

      aa_size = 1
      if m = alt.match(/(Indel|FrameShift)\((.*)\)/i)
        type = m[1]
        aas = m[2]
        aas_a = aas.split("")
        while aas_a.first == '-'
          mut_sequence[pos.to_i - 1] = ""
          aas_a.shift
        end
        if type == "Indel"
          aa_size = aas_a.length
          mut_sequence[pos.to_i - 1] = (aas_a * "") + mut_sequence[pos.to_i - 1]  if aas_a
        else
          mut_sequence[pos.to_i - 1..-1] = aas_a * ""
        end
      else
        type = 'SNV'
        mut_sequence[pos.to_i - 1] = alt
      end

      sizes.each do |size|
        size = size.to_i
        extra = (aa_size - 1)
        (size + extra).times do |i|
          start = pos.to_i - 1 - i + (aa_size - 1)
          next if start < 0

          eend = start + size - 1
          mut_epitope = mut_sequence[start..eend]
          mut_epitope = mut_epitope[0..-2] if mut_epitope[-1] == "*"
          next if mut_epitope.length < size
          next if mut_epitope.include?"*"

          wt_epitope1 = wt_sequence[start..eend]
          wt_epitope2 = wt_sequence[(start-extra)..(eend-extra)]
          wt_epitope = (i + 1) > ((size + extra)/2).round ?  wt_epitope1 : wt_epitope2
          wt_epitope = wt_epitope2 if wt_epitope.nil?

          wt_epitope = wt_epitope[0..-2] if wt_epitope && wt_epitope[-1] == "*"
          wt_epitope = nil if wt_epitope && wt_epitope.length < size
          wt_epitope = nil if wt_epitope && wt_epitope.include?("*")

          next if wt_epitope == mut_epitope
          epitopes << [wt_epitope, mut_epitope, i]
        end

        if type == "FrameShift"
          remain = mut_sequence.length - pos.to_i - size + 1
          remain.times do |i|
            start = pos.to_i + i
            eend = start + size - 1
            mut_epitope = mut_sequence[start..eend]
            mut_epitope = mut_epitope[0..-2] if mut_epitope[-1] == "*"
            next if mut_epitope.length < size
            next if mut_epitope.include?("*")
            epitopes << [nil, mut_epitope, start + 1]
          end
        end
      end


      next if epitopes.empty?
      [mi, Misc.zip_fields(epitopes)]
    end
  end

  dep :epitopes
  dep :epitopes_old
  task :cmp_epi => :array do
    tsv1 = step(:epitopes).load
    tsv2 = step(:epitopes_old).load
    diff = []
    tsv1.each do |k,v1|
      v2 = tsv2[k]
      p1 = v1[0].zip(v1[1]).sort_by{|a,b| b}
      p2 = v2[0].zip(v2[1]).sort_by{|a,b| b}
      if p1 != p2
        eee "diff" 
        wwww k
        wwww p1
        wwww p2
        wwww p1 - p2
        wwww p2 - p1
        diff << k
      end
    end
    diff
  end

end
