module Sample

  dep :genomic_mutations
  dep :organism
  dep Sequence, :reference, :positions => :genomic_mutations, :organism => :organism
  dep :genomic_mutation_annotations
  dep :sequence_ontology
  dep :DbNSFP, :principal => true
  dep :DbNSFP_pred, :principal => true
  dep :neo_epitopes, :principal => true
  dep Sequence, :genes, :positions => :genomic_mutations, :organism => :organism
  dep :expanded_vcf
  task :rtoledo_maf => :tsv do
    fields =<<-EOF.split("\n")
dbSNP
Ref
Alt
Type
Filter
Gene
Transcript
Exon
DNA_change
protein_change
CGI_input
SIFT_score
SIFT_pred
Polyphen2_HDIV_score
Polyphen2_HDIV_pred
Polyphen2_HVAR_score
Polyphen2_HVAR_pred
FATHMM_score
FATHMM_pred
gnomAD_exomes_AF
gnomAD_exomes_AC
gnomAD_exomes_AN
clinvar_clnsig
clinvar_trait
Cancer_Gene_Census
Biocarta_Pathway
KEGG_Pathway
MAF_tDNA
MAF_cfDNA
MAF_primary
    EOF

    organism = self.recursive_inputs[:organism]
    mutations = step(:genomic_mutations).load
    tsv = TSV.setup(mutations, :key_field => "Genomic Mutation", :fields => [], :type => :double, :namespace => organism)
    tsv.identifiers = Organism.identifiers(organism)

    annotations = step(:genomic_mutation_annotations).load
    consequence = step(:genomic_mutation_consequence).load
    reference = step(:reference).load
    reference.key_field = tsv.key_field

    tsv = tsv.attach annotations, :fields => ["DbSNP:RS ID", "GnomAD:AF"]
    tsv = tsv.attach reference, :fields => ["Reference Allele"]
    tsv.add_field "Alternative Allele" do |mutation|
      mutation.split(":")[2]
    end

    tsv = tsv.attach step(:sequence_ontology), :fields => ["SO Term"]

    genes = step(:genes).load.swap_id("Ensembl Gene ID", "Associated Gene Name")
    coding = Organism.gene_biotype(organism).tsv.change_key("Associated Gene Name").select("Biotype" => "protein_coding").keys
    genes.key_field = tsv.key_field

    tsv = tsv.attach genes, :fields => ["Associated Gene Name"]
    tsv.process "Associated Gene Name" do |genes|
      next if genes.nil?
      good = genes & coding
      good.any? ? good : genes
    end

    tsv.add_field "Mutated Isoform" do |mutation|
      mis = consequence[mutation] || []
    end

    tsv.add_field "Ensembl Protein ID" do |mutation|
      mis = consequence[mutation] || []
      mis.collect{|mi| mi.split(":").first}
    end

    tsv.add_field "AA Change" do |mutation|
      mis = consequence[mutation] || []
      mis.collect{|mi| mi.split(":").last}
    end

    tsv.attach step(:DbNSFP_pred).load, :fields => %w(SIFT_pred Polyphen2_HDIV_pred Polyphen2_HVAR_pred FATHMM_pred MetaSVM_pred)

    neo_epi = step(:neo_epitopes).load
    neo_epi_fixed = neo_epi.annotate({})
    neo_epi.each do |m,v|
      neo_epi_fixed[m.sub(/^chr/,'')] = v
    end
    tsv.attach neo_epi_fixed, :fields => ["MHCflurry MT Score", "MHCflurry WT Score"]


    tsv = Rbbt.with_workflow "Enrichment" do
      %w(kegg go_bp go_mf).each do |database| 
        enrichment = Enrichment.job(:enrichment, nil, :database => database, :list => genes.values.flatten, :organism => organism, :cutoff => 1.1, :fdr => false, :min_support => 0).run
        Log.tsv enrichment
        db_field = enrichment.key_field
        tsv.attach enrichment, :fields => [db_field]
        if tsv.fields.include? "GO ID"
          tsv.process "GO ID" do |go|
            GO.id2name go
          end 
          tsv.fields = tsv.fields.collect{|f| f == "GO ID" ? "GO ID (#{ database })" : f}
        end

      end
      tsv.swap_id "KEGG Pathway ID", "Pathway Name", :identifiers => KEGG.pathways
    end

    subset = case 
             when self.clean_name.include?('p136')
               "tumor"
             when self.clean_name.include?('IMN')
               "cfDNA"
             else
               "AltBaj"
             end

    exp_vcf = step(:expanded_vcf).load
    exp_vcf.fields = exp_vcf.fields.collect{|f| f.include?("normal:") ? f.sub("normal:", "germline_#{subset}:") : f }
    freq_fields = exp_vcf.fields.select{|f| %w(AD AF DP).include? f.split(":").last}
    tsv.attach exp_vcf, :fields => freq_fields


    tsv
  end
end
