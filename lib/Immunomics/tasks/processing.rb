require 'rbbt/sources/organism'
module Immunomics

  def self.netchop(file)
    cuts = {}
    TSV.traverse CMD.cmd("#{Rbbt.software.opt["netchop-3.1"].netchop.find} '#{file}'", :pipe => true), :type => :line do |line|
      next unless (line.include?("wt|") || line.include?("mu|")) && ! line.include?("Number of")
      parts = line.split(/\s+/)
      next unless parts[3] == "S"
      digest = parts[5][0..-1]
      cuts[digest] ||= []
      cuts[digest] << parts[1]
    end
    cuts
  end

  def self.tapmat_pred_fsa(file, size)
    scores = {}
    TSV.traverse CMD.cmd("#{Rbbt.software.opt["netCTLpan-1.1"].Linux_x86_64.bin.tapmat_pred_fsa.find} -mat '#{Rbbt.software.opt["netCTLpan-1.1"].data["tap.logodds.mat"].find}' '#{file}' -l #{size}", :pipe => true), :type => :line do |line|
      next if line =~ /^#/
      parts = line.split(/\s+/)
      scores[parts[2]] = parts[3]
    end
    scores
  end

  dep :mutation_flanking_sequence, :compute => :produce
  task :netchop => :tsv do 

    translations = {}
    io = TSV.traverse step(:mutation_flanking_sequence), :into => :stream  do |mi,values|
      lost, pos, n_flank, change, c_flank = values
      wt_sequence = [n_flank, lost, c_flank] * ""
      mut_sequence = [n_flank, change, c_flank] * ""

      md5 = Misc.digest(mi)[0..5]
      while translations.include?(md5)
        md5 = Misc.digest(mi + md5)[0..5]
      end
      translations[md5] = [mi, pos, n_flank.length]
      <<-EOF
>wt|#{md5}|mi|#{mi}
#{wt_sequence}
>mu|#{md5}|mi|#{mi}
#{mut_sequence}
      EOF
    end

    Open.write file('fasta'), io
    
    cuts = Immunomics.netchop(file('fasta'))

    tsv = TSV.setup({}, "Mutated Isoform~WT Cut position,Mut Cut positions#:type=:double")
    cuts.each do |digest,cs|
      type, md5 = digest.split("|")
      mi, pos, flank_size = translations[md5]
      cpos = cs.collect{|c| c.to_i + pos.to_i - 1 - flank_size} 
      tsv[mi] ||= [[],[]]
      if type == "wt"
        tsv[mi][0] = cpos
      else
        tsv[mi][1] = cpos
      end
    end
    tsv
  end

  dep :mutation_flanking_sequence, :compute => :produce
  input :tap_sizes, :array, "Peptide size", (8..16).to_a
  task :tapmat_pred_fsa => :tsv do |sizes|

    translations = {}
    io = TSV.traverse step(:mutation_flanking_sequence), :into => :stream  do |mi,values|
      lost, pos, n_flank, change, c_flank = values
      wt_sequence = [n_flank, lost, c_flank] * ""
      mut_sequence = [n_flank, change, c_flank] * ""

      wt_sequence.gsub!('-', '')
      mut_sequence.gsub!('-', '')

      wt_sequence.gsub!('*', '')
      mut_sequence.gsub!('*', '')

      md5 = Misc.digest(mi)[0..5]
      while translations.include?(md5)
        md5 = Misc.digest(mi + md5)[0..5]
      end

      translations[md5] = [mi, pos, n_flank.length]
      <<-EOF
>wt|#{md5}|mi|#{mi}
#{wt_sequence}
>mu|#{md5}|mi|#{mi}
#{mut_sequence}
      EOF
    end

    Open.write file('fasta'), io
    
    scores = TSV.setup({}, "Peptide~Score#:type=:single#:cast=:to_f")
    sizes.each do |size|
      scores.merge!(Immunomics.tapmat_pred_fsa(file('fasta'), size))
    end

    scores
  end

  #helper :best_tap_old do |matrix,epi,cuts=nil,cpos=nil|
  #  scores = {}
  #  matrix.each do |fragment,score|
  #    if fragment.end_with?(epi)
  #      l = fragment.length
  #      next if cuts && cpos && ! cuts.include?((cpos - l).to_s)
  #      scores[l] ||= []
  #      scores[l] << score
  #    end
  #  end
  #  scores.values.flatten.max
  #end

  #helper :best_tap do |matrix,epi,cuts=nil,cpos=nil|
  #  scores = {}
  #  cuts = Set.new cuts.collect{|c| c.to_i} if cuts
  #  max = nil
  #  matrix.each do |fragment,score|
  #    next if max && score <= max
  #    if fragment.end_with?(epi)
  #      l = fragment.length
  #      next if cuts && cpos && ! cuts.include?(cpos - l)
  #      max = score
  #    end
  #  end
  #  max
  #end

  helper :best_tap_index do |index,epi,cuts=nil,cpos=nil|
    scores = {}
    cuts = Set.new cuts.collect{|c| c.to_i} if cuts
    max = nil
    matches = index.prefix(epi.reverse)
    matches.each do |fragment|
      score = index[fragment]
      next if max && score <= max
      next if cuts && cpos && ! cuts.include?(cpos - fragment.length)
      max = score
    end
    max
  end

  dep :epitopes
  dep :netchop
  dep :tapmat_pred_fsa
  task :processed_epitopes => :tsv do
    netchop = step(:netchop).load
    tapmat_pred_fsa = step(:tapmat_pred_fsa).load

    inverse_index = file(:index)

    tapmat_pred_fsa_index = Persist.open_tokyocabinet(inverse_index, true, :float, TokyoCabinet::BDB)

    tapmat_pred_fsa.with_monitor do
      tapmat_pred_fsa_index.write
      tapmat_pred_fsa.through do |p,s|
        tapmat_pred_fsa_index[p.reverse] = s
      end
    end

    tapmat_pred_fsa_index.read

    netchop.unnamed = true
    tapmat_pred_fsa.unnamed = true

    parser = TSV::Parser.new step(:epitopes)
    fields = parser.fields
    dumper = TSV::Dumper.new parser.options.merge(:fields => fields + ["WT Cut", "Mut Cut", "WT Truncated", "Mut Truncated", "WT TAP", "Mut TAP", "WT TAP Strict", "Mut TAP Strict"])
    dumper.init
    TSV.traverse parser, :into => dumper, :bar => true, :cpus => 10 do |mi, values|
      mi = mi.first if Array === mi
      prot, change = mi.split(":")
      wtaa, pos, altaa = change.partition /\d+/

      zipped_res = Misc.zip_fields(values).collect do |wt,mut,off|
        if ! wt.empty?
          wtcpos = pos.to_i + off.to_i
        else
          wtcpos = nil
        end

        if ! mut.empty?
          mutcpos = pos.to_i + off.to_i
        else
          mutcpos = nil
        end

        wt_cuts, mut_cuts = netchop[mi]
        wt_cuts = [] if wt_cuts.nil?
        mut_cuts = [] if mut_cuts.nil?

        wt_good = wt_cuts.include?(wtcpos.to_s)
        mut_good = mut_cuts.include?(mutcpos.to_s)

        wt_positions = ((wtcpos - wt.length + 1)..(wtcpos - 1)).to_a.collect{|p| p.to_s} if wtcpos
        mut_positions = ((mutcpos - mut.length + 1)..(mutcpos - 1)).to_a.collect{|p| p.to_s} if mutcpos

        wt_truncated = (wt_cuts & wt_positions).any? if wtcpos
        mut_truncated = (mut_cuts & mut_positions).any? if mutcpos

        wt_tap = best_tap_index tapmat_pred_fsa_index, wt unless wt.empty?
        mut_tap = best_tap_index tapmat_pred_fsa_index, mut unless mut.empty?

        wt_tap_strict = best_tap_index tapmat_pred_fsa_index, wt, wt_cuts, wtcpos unless wt.empty?
        mut_tap_strict = best_tap_index tapmat_pred_fsa_index, mut, mut_cuts, mutcpos unless mut.empty?

        [wt_good, mut_good, wt_truncated, mut_truncated, wt_tap, mut_tap, wt_tap_strict, mut_tap_strict]
      end

      [mi, values + Misc.zip_fields(zipped_res)]
    end
  end

end

