module Sample

  %w(epitopes netchop tapmat_pred_fsa processed_epitopes mhcFlurry).each do |task_name|
    dep :organism
    dep :mi, :organism => :organism
    dep_task task_name.to_sym, Immunomics, task_name.to_sym, :mutated_isoforms => :mi, :organism => :organism
  end

  dep :processed_epitopes, :compute => :produce
  dep :mhcFlurry, :compute => :produce
  task :epitope_info => :tsv do
    processed_epitopes = step(:processed_epitopes).load
    mhcFlurry_mut = step(:mhcFlurry).path.tsv :key_field => "Mutated epitope", :merge => true, :type => :double
    mhcFlurry_wt = step(:mhcFlurry).path.tsv :key_field => "Wildtype epitope", :merge => true, :type => :double

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

    processed_epitopes
  end
end
