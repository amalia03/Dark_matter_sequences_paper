#!/bin/bash

echo "Linear regression comparisons of FASTA kmers" > regression_results.txt

aligned="/home/ama/data/miseq_reassembly/orf_analysis/miseq_samples/merge_read_processing/all_aligned_300.fa"
unaligned="/home/ama/data/miseq_reassembly/orf_analysis/miseq_samples/merge_read_processing/all_unaligned_300.fa"
cap_C="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/kmer_analysis/group_comparison/C1C_worm.fa"
scol_C="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/kmer_analysis/group_comparison/S1C_worm.fa"
geans120A_head="/DATA/MarineGenomics/ama/data/GEANS/120A_head_nless.fa"
geans120A_tail="/DATA/MarineGenomics/ama/data/GEANS/120A_tail_nless.fa"
geans120B_head="/DATA/MarineGenomics/ama/data/GEANS/120B_head_nless.fa"
geans840A_head="/DATA/MarineGenomics/ama/data/GEANS/840A_head_nless.fa"

human_dna="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/kmer_analysis/group_comparison/human_seqs/Homo_sapiens_dna_sansn.fa"
human_cds="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/kmer_analysis/group_comparison/human_seqs/Homo_sapiens_cds_sansn.fa"

zebraf_dna="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/kmer_analysis/group_comparison/danio_rerio/Danio_rerio.GRCz11.dna.chromosome.sansn.fa"
zebraf_cds="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/kmer_analysis/group_comparison/danio_rerio/Danio_rerio.GRCz11.cds.sans.fa"
lumpfish="/home/ama/data/blastdb/lumpfish_assembly/lumpfish_trin_assembly.fa"

file_array=($cap_C \
    $scol_C \
    $geans120A_head \
    $geans120A_tail \
    $geans120B_head \
    $geans840A_head \
    $human_dna \
    $human_cds \
    $zebraf_dna \
    $zebraf_cds \
    $lumpfish
)


#For aligned vs unaligned

echo "Aligned Vs Unaligned" >> regression_results.txt
echo "Aligned Vs Unaligned"
./kmer_regression_function.R $aligned $unaligned >> regression_results.txt

for i in ${file_array[@]}
do

id="$(ls $i | rev | cut -d "/" -f 1 | rev )"
##For Aligned
echo "Aligned vs $id" >> regression_results.txt
echo "Aligned vs $id"
./kmer_regression_function.R $aligned $i >> regression_results.txt


##For Unaligned
echo "Unaligned vs $id" >> regression_results.txt
echo "Unaligned vs $id"
./kmer_regression_function.R $unaligned $i >> regression_results.txt

done
