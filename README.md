# Dark_matter_sequences_section
In this repository you will find scripts that were used for the research paper "The dark matter of environmental transcriptomics: what information can we obtain from priorly unannotated benthic mRNA sequences?" (currently unpublished). Below is the process we follow to analyse our data. Some scipts are located in different depositories and noted as such. 

## 1. Get longest sequences, BLAST locat alignment to NCBI and partition the sequences to mapped and unmapped. 

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

See README.md file in kmer analysis folder.

## 3. ORF analysis

For the Open reading frame analysis we created two simulated dataset types to compare the ORFs against our aligned and unaligned reads. Then we used the Orf prefictor program to identify orfs in our sequences and then a cusotm R script to analyse the results. 

For simulating reads, we first downloaded genomes from model organisms. The zebrafish (Danio rerio) genome will be used as an example here. 
The genome was downloaded from Ensembl using the command 

`Dark_matter_sequences_paper/orf_section/wget_danio.sh`


## 4. Protein discovery analysis

