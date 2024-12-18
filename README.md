# Dark_matter_sequences_section
In this repository you will find scripts that were used for the research paper "The dark matter of environmental transcriptomics: what information can we obtain from priorly unannotated benthic mRNA sequences?" (currently unpublished). Below is the process we follow to analyse our data. Some scipts are located in different depositories and noted as such. 

## 1. Get longest sequences, BLAST locat alignment to NCBI and partition the sequences to mapped and unmapped. 

1. Sequences were merged using PEAR (`make_match_unmatch_files/pear_pairs.sh`).
2. Longer merged sequences were kept (>= 300 bp) (`make_match_unmatch_files/filter_fasta_lengths.R`).
3. Merged sequences were BLASTed against the complete NCBI database using very low thresholds (no pid, evalue e-5).
4. The total sequences were separated in two FA files, one with sequences that had matches to BLAST and another that had none(`make_match_unmatch_files/fetch_unaligned_seqs.pl`, `make_match_unmatch_files/fetch_aligned_seqs.pl`)
   
This last two output files were the ones that were used for the rest of the analysis below. 

## 2. Kmer analysis

See README.md file in kmer analysis folder.

## 3. ORF analysis

For the Open reading frame analysis we created two simulated dataset types to compare the ORFs against our aligned and unaligned reads. Then we used the Orf prefictor program to identify orfs in our sequences and then a cusotm R script to analyse the results. 
Orfpredictor: Min, X.J., Butler, G., Storms, R. and Tsang, A. OrfPredictor: predicting protein-coding regions in EST-derived sequences. Nucleic Acids Res., 2005

-----
For simulating dataset based on model genomes: 
1. we first downloaded genomes from model organisms. The zebrafish (Danio rerio) genome will be used as an example here. 
The genome was downloaded from Ensembl using the command `orf_section/wget_danio.sh`.

2. The genomes were then concatenated together, while also removing all N segments using `orf_section/cat_nless_genome.pl`

3. The large sequence file is then cut into random segments using the same lengths found in the original merged files (`orf_section/rando_genome_snipper.sh`). 

4. Finally, OrfPredictor was used through the following script
`perl /home/ama/bin/OrfPredictor/OrfPredictor.pl zebrafish_sim_300bp.fa blank.txt 0 both 1 300b`

------

For the simulated dataset based on random reads:

1. Determine the total nucleotide percent of the merged sequence fasta file using the command `orf_section/find_nucl_percent.pl`:

2. Create random sequences based on the length of the merged sequence files as well as the expected nucleotide proportions provided by the previous step `orf_section/find_nucl_percent.pl`.

3. Then I used a loop version of the OrfPredictor `orf_section/orfpredictor_command.sh`.

----

Orf regions for aligned and unaligned sets were found using Orfpredictor, and outputs from the generated genome orfs, the randomised dataset orfs and the aligned and unaligned orfs were anlaysed and plotted using `orf_section/codon_usage.R`.

## 4. Protein discovery analysis

For this analysis, we prepared 5 datasets in addition to the original BLAST search against the NCBI database, that is: 

1. A BLAST search against the nr database.
2. A BLASTx search against the Uniprot database
3. A search using "Infernal" against Rfam to search for non-coding rna sequences.
4. An Interpro search against the InterPro database
5. An MG-RAST run

-----

### Rfam dataset preparation

1. Downladed the RFAM database.
2. Run the the cmpress command (`db_comparison/cmpress.sh` ) from Infernal to format the downloaded cm files
3. Run the cmscan command on the cm files.

For the analysis later I only use a tabulated array of all the ids that matched to RFAM. 

----

### Interpro datasets preparations 

We first used the ORFs from the previous section, used aligned and unaligned separately. 

1. We split the sequences into files using the `/db_comparison/select_orfs_by_frame.pl` command.

And executing it using the bash command:

`./select_orfs_by_frame.pl all_unaligned_300.fa unaligned_orfs/ORF6frame.txt 0 orf_pos orf_neg 500 > orf_stats.tsv`
`./select_orfs_by_frame.pl all_aligned_300.fa unaligned_orfs/ORF6frame.txt 0 orf_pos orf_neg 500 > orf_stats.tsv`

where $nuc_file is the FASTA file, $pep_file is the peptide file $min_length is the minimum length,  $pos_odir is the positive strand directory, $neg_odir is the negative strand directory and , $group_size is the number of sequences per file. 

The output creates fasta files that in this case contain 500 fasta sequences from longest to shortest sequences.

2. We ran InteProScan v5. (interproscan.sh in bin directory) into the two strand directories using the `db_comparison/run_scans.sh`
command. 

Which activates the `db_comparison/run_interproscan.sh` command that has been included in the two strand directories (probably there's a simpler way of doing it). 

------

A protein database comparison was then made using `db_comparison/db_comparison_analysis.R`  script. This script outputs an Upsett plot of db matches as well a grid that shows how much exlusive and shared information each database provides relative to all others

`db_comparison/db_comparison_analysis.R`

### Interpro Analysis

Data analysis and vizualizasions on the InterPro database results were made using the following `db_comparison/interpro_analysis.R` R script. 

Which should be paired with the `db_comparison/pep_functions.R` that contain some major functions uses in the main script. 

This is the most detailed script of the bunch. It is used to generate three major outputs that were used for the paper. 
The 3 main plots generated are i) two barplots showing alignment results to the InterPro databases, ii) a protein scatterplot indicated the abundance of the different protein groups that aligned within each database and finally ii) an upsett plot indicating the matches between each group. 

It also includes some initial venn diagram plots that I used to compare the results between the databases, and a tabular output of the highest proteins and their descriptions. 
