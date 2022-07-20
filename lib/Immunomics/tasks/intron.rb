module Immunomics
  dep Sequence, :genes
  task :intron_mutations => :tsv do
    organism = recursive_inputs[:organism]
    TSV.traverse step(:genes) do |mutation,gene|
      chr
    end
  end
end
