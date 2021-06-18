module Sample

  dep :stringtie
  input :abundancy_field, :select, "Abundancy field to extract from stringtie", "FPKM", :select_options => %w(FPKM TMP)
  task :gene_abundancies => :tsv do |abundancy_field|

    dumper = TSV::Dumper.new :key_field => "Ensembl Transcript ID", :fields => [abundancy_field], :type => :single, :cast => :to_f 
    dumper.init
    TSV.traverse step(:stringtie), :into => dumper, :type => :array do |line|
      next if line =~ /^#/
      parts = line.split("\t")
      info = parts[8]

      gene = transcript = reference = abundance = nil
      info.split(/;\s+/).each do |i|
        key, value = i.split(" ")
        case key
        when "gene_id"
          gene = value.gsub('"', '')
        when "transcript_id"
          transcript = value.gsub('"', '')
        when "reference_id"
          reference = value.gsub('"', '')
        when abundancy_field
          abundance = value.gsub('"', '').to_f
        end
      end

      next if abundance.nil? || reference.nil?
      [reference, abundance]
    end
  end

  dep :gene_abundancies
  dep :organism
  dep :mi, :organism => :organism
  dep_task :mi_expression, Immunomics, :mi_expression, :gene_abundancies => :gene_abundancies, :organism => :organism, :mutated_isoforms => :mi

  dep :genomic_mutations
  dep :RNA_BAM
  dep_task :genomic_mutation_rna_evidence, HTS, :mutation_support, :mutations => :genomic_mutations, :bam => :RNA_BAM

  dep :genomic_mutation_rna_evidence, :compute => :produce
  dep :genomic_mutation_consequence, :compute => :produce
  task :mi_rna_evidence => :tsv do

    evidence = step(:genomic_mutation_rna_evidence).load
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => evidence.fields, :type => evidence.type
    dumper.init
    TSV.traverse step(:genomic_mutation_consequence), :into => dumper do |mut,mis|
      mut = mut.first if Array === mut
      value = evidence[mut]
      next if value.nil?
      res = mis.collect{|mi| [mi, value] }

      res.extend MultipleResult

      res
    end

  end

  dep :pyclone, :compute => :produce, :canfail => true
  dep :genomic_mutation_consequence, :compute => :produce
  task :mi_prevalence => :tsv do
    pyclone = step(:pyclone).load
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["cellular_prevalence"], :type => :single, :cast => :to_f
    dumper.init
    TSV.traverse step(:genomic_mutation_consequence), :into => dumper do |mut,mis|
      mut = mut.first if Array === mut
      value = pyclone["b'#{mut}'"]
      next if value.nil?
      value = value["cellular_prevalence"]
      res = mis.collect{|mi| [mi, value] }

      res.extend MultipleResult

      res
    end

  end

  %w(epitopes netchop tapmat_pred_fsa processed_epitopes mhcFlurry).each do |task_name|
    dep :organism
    dep :mi, :organism => :organism
    dep_task task_name.to_sym, Immunomics, task_name.to_sym, :mutated_isoforms => :mi, :organism => :organism
  end

  dep :processed_epitopes, :compute => :produce
  dep :mhcFlurry, :compute => :produce
  dep :mi_expression, :compute => :produce, :canfail => true
  dep :mi_rna_evidence, :compute => :produce, :canfail => true
  dep :mi_prevalence, :compute => :produce, :canfail => true
  task :epitope_info => :tsv do
    processed_epitopes = step(:processed_epitopes).load
    mhcFlurry_mut = step(:mhcFlurry).path.tsv :key_field => "Mutated epitope", :merge => true, :type => :double
    mhcFlurry_wt = step(:mhcFlurry).path.tsv :key_field => "Wildtype epitope", :merge => true, :type => :double

    mi_expression = step(:mi_expression).load if  step(:mi_expression).done?
    mi_rna_evidence = step(:mi_rna_evidence).load if  step(:mi_rna_evidence).done?

    pyclone = step(:mi_prevalence).load if  step(:mi_prevalence).done?

    fields = mhcFlurry_mut.fields.select{|f| f.include? "mhcflurry" }
    mut_fields = fields.select{|f| f.include? "Mutated" }
    wt_fields = fields - mut_fields

    mut_fields.each do |field|
      processed_epitopes.add_field field do |k,values|
        mutated = values["Mutated"] || []
        mutated.collect{|e| mhcFlurry_mut.include?(e)? Misc.min(mhcFlurry_mut[e][field]) : nil }
      end
    end

    wt_fields.each do |field|
      processed_epitopes.add_field field do |k,values|
        wt = values["Wildtype"]
        wt.collect{|e| mhcFlurry_wt.include?(e)? Misc.min(mhcFlurry_wt[e][field]) : nil }
      end
    end

    if mi_expression
      processed_epitopes.add_field mi_expression.fields.first do |k,values|
        [mi_expression[k] || 0] * values.first.length
      end
    end

    if mi_rna_evidence
      processed_epitopes.add_field mi_rna_evidence.fields.first do |k,values|
        [mi_rna_evidence[k].nil? ? 0 : mi_rna_evidence[k].first] * values.first.length
      end
    end

    if mi_rna_evidence
      processed_epitopes.add_field mi_rna_evidence.fields.last do |k,values|
        [mi_rna_evidence[k].nil? ? 0 : mi_rna_evidence[k].last] * values.first.length
      end
    end

    if pyclone
      processed_epitopes.add_field "cellular_prevalence" do |k,values|
        [pyclone[k].nil? ? 0 : pyclone[k]] * values.first.length
      end
    end




    processed_epitopes
  end


end
