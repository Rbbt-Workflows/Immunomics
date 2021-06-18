module Immunomics

  input :mutated_isoforms, :array, "Mutated Isoforms"
  input :gene_abundancies, :tsv, "Gene (or transcript) abundancies"
  input :match_expression, :select, "Match epitopes by protein isoform and transcript or by gene", 'gene', :select_options => %w(gene transcript)
  input :organism, :string, "Organism code", Organism.default_organism("Hsa")
  task :mi_expression => :tsv do |mis,gene_abundancies,match_expression,organism|
    ensp2enst = Organism.transcripts(organism).tsv :key_field => "Ensembl Protein ID", :fields => "Ensembl Transcript ID", :type => :single, :persist => true
    enst2ensg = Organism.transcripts(organism).tsv :key_field => "Ensembl Transcript ID", :fields => "Ensembl Gene ID", :type => :single, :persist => true

    Log.tsv gene_abundancies
    gene_abundancies = gene_abundancies.change_key "Ensembl Gene ID", :identifiers => Organism.transcripts(organism)  if match_expression == 'gene'

    dumper = TSV::Dumper.new :key_field => "Epitope", :fields => [gene_abundancies.fields.first], :type => :flat, :cast => :to_f
    dumper.init
    TSV.traverse mis, :into => dumper do |mi|
      protein = mi.split(":").first
      transcript = ensp2enst[protein]

      if match_expression == 'gene'
        gene = enst2ensg[transcript]
        abundance = gene_abundancies[gene]
      else
        abundance = gene_abundancies[transcript]
      end

      next if abundance.nil?

      [mi, [abundance]]
    end

    tsv = TSV.open dumper.stream, :merge => true
    tsv.to_single{|l| Misc.sum(l) }
  end

end
