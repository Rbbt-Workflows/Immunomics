require 'rbbt/sources/organism'
module Immunomics

  CMD.tool 'mhcflurry-predict' do
    CMD.cmd_log('pip install mhcflurry && mhcflurry-downloads fetch models_class1_presentation')
  end
  
  dep :epitopes
  input :alleles, :array, "List of alleles", nil, :required => true
  task :mhcFlurry => :tsv do |alleles|
    str = 'allele,peptide' << "\n"
    TSV.traverse step(:epitopes), :into => str do |mi,values|
      mi = mi.first if Array === mi
      wt_epitopes, mut_epitopes, offsets = values
      next if mut_epitopes.nil? or mut_epitopes.compact.empty?
      all_epitopes = (wt_epitopes + mut_epitopes).compact.uniq
      all_epitopes.reject!{|e| e.empty?}

      res = []

      res = alleles.collect{|a| all_epitopes.compact.uniq.collect{|epitope| [a,epitope] * "," } }.flatten

      res * "\n" + "\n"
    end

    input = file('input.csv')
    output = file('output.csv')

    Open.write(input, str)
    CMD.cmd_log('mhcflurry-predict'," #{ input } --out #{output}" )

    mhcflurry_scores = {}
    TSV.traverse output, :type => :array do |line|
      next if line =~ /mhcflurry_prediction/
      allele, epitope, pred, percent, processing, presentation = line.split(",")
      mhcflurry_scores[[allele, epitope]] = [pred, percent, processing, presentation]
    end

    dumper = TSV::Dumper.new :key_field => "Mutation ID", :type => :list, :fields => ["Mutated Isoform", "Allele", "Mutated epitope", "Wildtype epitope", "Offset",
                                                                      "Mutated mhcflurry prediction", "Mutated mhcflurry percentile", "Mutated mhcflurry processing", "Mutated mhcflurry presentation", 
                                                                      "Wildtype mhcflurry prediction", "Wildtype mhcflurry percentile", "Wildtype mhcflurry processing", "Wildtype mhcflurry presentation"]
    dumper.init
    TSV.traverse step(:epitopes), :into => dumper do |mi,values|
      mi = mi.first if Array === mi
      wt_epitopes, mut_epitopes, offsets = values
      next if mut_epitopes.nil? or mut_epitopes.compact.empty?

      res = []
      Misc.zip_fields(values).each do |wt_epitope,mut_epitope,offset|
        alleles.each do |allele|
          size = mut_epitope.length
          id = [mi, size, offset, allele] * "."
          wt_scores = mhcflurry_scores[[allele, wt_epitope]]
          mut_scores = mhcflurry_scores[[allele, mut_epitope]]
          wt_scores = [nil] * mut_scores.length if wt_scores.nil?
          res << [id, [mi, allele, mut_epitope, wt_epitope, offset] + mut_scores + wt_scores]
        end
      end
      res.extend MultipleResult
      res 
    end
  end

end
