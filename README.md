# Dark_matter_sequences_section
In this repository you will find scripts that were used for the research paper "The dark matter of environmental transcriptomics: what information can we obtain from priorly unannotated benthic mRNA sequences?" (currently unpublished). Below is the process we follow to analyse our data. Some scipts are located in different depositories and noted as such. 

## 1. Keeping the longest sequences and partitioning the sequences to mapped and unmapped. 

1. Sequences were merged using PEAR.
   `Dark_matter_sequences_paper/make_match_unmatch_files/pear_pairs.sh`
3. Longer merged sequences were kept (>= 300 bp).
   `Dark_matter_sequences_paper/make_match_unmatch_files/filter_fasta_lengths.R`
5. Merged sequences were BLASTed against the complete NCBI database using very low thresholds (no pid, evalue e-5).
6. The total sequences were separated in two FA files, one with sequences that had matches to BLAST and another that had none.
   `Dark_matter_sequences_paper/make_match_unmatch_files/fetch_unaligned_seqs.pl`
   `Dark_matter_sequences_paper/make_match_unmatch_files/fetch_aligned_seqs.pl`
   
This last two output files were the ones that were used for the rest of the analysis below. 

## 2. Kmer analysis
## 3. ORF analysis
## 4. Protein discovery analysis

