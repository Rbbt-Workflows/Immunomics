module Immunomics

  dep Sample, :mi, :vcf => true
  dep_task :vcf_epitopes, Immunomics, :epitopes, :mutated_isoforms => :mi

end
