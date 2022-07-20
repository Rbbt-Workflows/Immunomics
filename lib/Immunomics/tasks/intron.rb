module Immunomics

  input :mutations, :array, "Mutations"
  input :organism, :string, "Organism", Organism.default_code('Hsa')
  input :flank_size, :integer, "Flank size", 100
  task :mutation_dna_sequences => :tsv do |mutations,organism,flank_size|
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Pre", "Reference", "Post"], :type => :list, :namespace => organism
    dumper.init
    chromosome_files = {}

    TSV.traverse mutations, :bar => self.progress_bar("Mutation DNA sequences"), :type => :array, :into => dumper do |mutation|
      begin
        chr, pos, alt = mutation.split(/[\s:\t]+/)
        next if pos.nil?
        chr.sub!(/^chr/i,'')
        chr = "MT" if chr == "M"
        file = chromosome_files[chr] ||= begin
                                           Sequence.chromosome_file(organism, chr)
                                         rescue Exception
                                           :missing
                                         end
        next if file == :missing

        pos = pos.to_i

        file.seek pos - flank_size

        alt_length = 1
        alt_length = alt.scan("-").length if alt.include?("-")

        pre = file.read(flank_size)
        ref = file.read(alt_length)
        post = file.read(flank_size)

        [mutation, [pre, ref, post]]
      rescue Exception
        Log.exception $!
        [mutation, "?"]
      end
    end
  end

  helper :compare_aa_seqs do |wt_translations,mut_translations|
    res = []
    wt_translations.zip(mut_translations).each_with_index do |p,i|
      wt, mut = p
      wt = wt.chars
      mut = mut.chars
      wtr = wt.reverse
      mutr = mut.reverse

      common_pre = 0
      wt.each_with_index do |c,i|
        break if mut[i] != c
        common_pre += 1
      end

      common_post = 0
      wtr.each_with_index do |c,i|
        break if mutr[i] != c
        common_post += 1
      end

      left = wt[0..common_pre-1] * ""
      right = wtr[0..common_post-1].reverse * ""
      lost = (wt[common_pre..-(common_post + 1)] || []) * ""
      change = (mut[common_pre..-(common_post + 1)] || []) * ""

      right, lost = lost, "" if common_post == 0

      values = [lost, common_pre + 1, left, change, right]

      res << values
    end

    res

  end

  input :sizes, :array, "Size of epitope", [9]
  dep :mutation_dna_sequences, :flank_size => :placeholder do |jobname,options|
    max_size = options[:sizes].collect{|s| s.to_i }.max
    {:inputs => options.merge(:flank_size => (max_size + 2) * 3)}
  end
  task :mutation_open_frames => :tsv do |sizes|
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation Reading Frame", :fields => ["Lost", "Position of Change", "N Flank", "Change", "C Flank"], :type => :list
    dumper.init
    TSV.traverse step(:mutation_dna_sequences), :into => dumper do |mutation,values|
      mutation = mutation.first if Array === mutation
      chr, pos, alt = mutation.split(":")
      pre, ref, post = values

      alt = alt.gsub("-", '')

      wildtype = pre + ref + post
      mutated = pre + alt + post

      wt_translations = []
      mut_translations = []
      wt_translations_rev = []
      mut_translations_rev = []
      3.times do |orf|
        wt_translations << Bio::Sequence::NA.new(wildtype[orf..-1]).translate
        mut_translations << Bio::Sequence::NA.new(mutated[orf..-1]).translate
        wt_translations_rev << Bio::Sequence::NA.new(wildtype[0..-(orf + 1)].reverse).complement.translate
        mut_translations_rev << Bio::Sequence::NA.new(mutated[0..-(orf + 1)].reverse).complement.translate
      end

      res = []
      res.extend MultipleResult

      compare_aa_seqs(wt_translations, mut_translations).each_with_index do |values,i|
        key = ["RF","+", i, mutation] * ":"
        next if values[0].empty? && values[3].empty?
        res << [key, values]
      end

      compare_aa_seqs(wt_translations_rev, mut_translations_rev).each_with_index do |values,i|
        key = ["RF","-", i, mutation] * ":"
        next if values.first.empty?
        res << [key, values]
      end

      res.extend MultipleResult
      res
    end
  end

  dep :mutation_open_frames
  dep_task :open_epitopes, Immunomics, :epitopes,
    "Immunomics#mutation_flanking_sequence" => :mutation_open_frames,
    :flank_size => nil, :mutated_isoforms => nil, :sizes => nil

  dep :open_epitopes
  dep Sequence, :genes, :positions => :mutations, :vcf => false
  dep Sequence, :mutated_isoforms_fast, :coding => true, :non_synonymous => false, :principal => false, :watson => true, :vcf => false
  task :open_intron_epitopes => :tsv do 
    Step.wait_for_jobs dependencies
    genes = step(:genes).load
    mis_muts = step(:mutated_isoforms_fast).load.keys

    organism = recursive_inputs[:organism]
    gene_strands = Organism.gene_positions(organism).tsv :fields => ["Strand"], :type => :single, :persist => true

    parser = TSV::Parser.new step(:open_epitopes)
    dumper = TSV::Dumper.new parser.options.merge(:fields => parser.fields + ["Ensembl Gene ID"])
    dumper.init
    TSV.traverse parser, :bar => self.progress_bar("Selecting intron mutations"), :into => dumper do |key,values|
      key = key.first if Array === key
      parts = key.split(":")
      strand = parts[1] == "+" ? "1" : "-1"
      mutation = parts.values_at(3,4,5) * ":"
      next if mis_muts.include? mutation
      next unless genes.include? mutation
      overlaps = genes[mutation]
      next unless gene_strands.values_at(*genes[mutation]).include? strand
      [key, values + [overlaps]]
    end
  end
end
