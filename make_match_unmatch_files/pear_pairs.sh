#!/bin/bash

#A script for initiating the default pear command for all our aligned pairs
read_loc="/DATA/MarineGenomics/lmj/genomes/ama/miseq_reassembly/trimmed_only"
#merged_loc="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/merging_reads_pear/"
##Edit, 7.1.22:changed the file location as the names have changed..
merged_loc="/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/process_reads_for_ORF/"

# pear -f 183AD_unal_r1.fq -r 183AD_unal_r2.fq -o 183AD
cd $read_loc
id="$(ls *fq.gz | sed 's/_trim_[1,2]P_tr.fq.gz//' | uniq)"
echo $id

r1="_trim_1P_tr.fq.gz"
r2="_trim_2P_tr.fq.gz"


for i in $id;
do

#echo "pear -f $unal_loc$i$r1 -r $unal_loc$i$r2 -o $merged_loc$i"

pear -f $unal_loc$i$r1 -r $unal_loc$i$r2 -o $merged_loc$i

done
