module Immunomics

  input :organism, :select, "Organism code", "Hsa/may2017", :select_options => {"Hsa/feb2014" => "hg19", "Hsa/may2017" => "hg38"}
  dep Sample, :mi, :vcf => true
  dep_task :vcf_epitopes, Immunomics, :epitopes, :mutated_isoforms => :mi

end
