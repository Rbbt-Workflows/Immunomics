module Immunomics

  input :vcf_file, :file, "VCF file"
  input :organism, :select, "Organism code", "Hsa/may2017", :select_options => {"Hsa/feb2014" => "hg19 (Hsa/feb2014)", "Hsa/may2017" => "hg38 (Hsa/may2017) "}
  dep Sequence, :mutated_isoforms_fast, :vcf => true, :mutations => :vcf_file, :organism => :organism, :watson => true, :coding => true, :non_synonymous => true
  task :vcf_mis => :array do |vcf,organism|
    step(:mutated_isoforms_fast).load.values.flatten.uniq
  end

  input :organism, :select, "Organism code", "Hsa/may2017", :select_options => {"Hsa/feb2014" => "hg19 (Hsa/feb2014)", "Hsa/may2017" => "hg38 (Hsa/may2017) "}
  dep :vcf_mis, :organism => :organism
  dep_task :vcf_epitopes, Immunomics, :epitopes, :mutated_isoforms => :vcf_mis, :flank_size => 20, :organism => :organism

  input :organism, :select, "Organism code", "Hsa/may2017", :select_options => {"Hsa/feb2014" => "hg19 (Hsa/feb2014)", "Hsa/may2017" => "hg38 (Hsa/may2017) "}
  dep Sequence, :genomic_mutations, :quality => nil, :filters => []
  dep_task :vcf_open_epitopes, Immunomics, :open_epitopes, :mutations => :genomic_mutations, :organism => :organism

  input :organism, :select, "Organism code", "Hsa/may2017", :select_options => {"Hsa/feb2014" => "hg19 (Hsa/feb2014)", "Hsa/may2017" => "hg38 (Hsa/may2017) "}
  dep Sequence, :genomic_mutations, :quality => nil, :filters => []
  dep_task :vcf_open_intron_epitopes, Immunomics, :open_intron_epitopes, :mutations => :genomic_mutations, :organism => :organism
end
