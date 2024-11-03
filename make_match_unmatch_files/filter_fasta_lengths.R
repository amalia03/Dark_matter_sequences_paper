#!/bin/R
## This script takes a fasta file as an input and only outputs sequences of a specific length. 
## Edit: The original version tabulated the sequences next to the id instead of below it, which meant I created an extra file and script to untabulate the files. 

#Import fasta
miseq.fa <- readLines("all_sequences_Illumina_Miseq.fasta")

#Seperate id and sequence lines (assuming the sequence is not formatted into multiple lines
seqs <- miseq.fa[seq(2,length(miseq.fa), 2)]
ids <- miseq.fa[seq(1,length(miseq.fa)-1, 2)]

##Make a variable with the length of the sequences
read.l <- nchar(seqs)
###
Processing time, in this case, keep sequences with length larger than 300bp.
seqs.new <- seqs[read.l>300]
ids.new <- ids[read.l>300]

#Make a new variable that in the fasta format
new.fa <- capture.output(for (i in 1:length(seqs.new)) {
                             cat(ids.new[i], "\n")  # Print each element from A
                             cat(seqs.new[i], "\n")  # Print each element from B
                         })
#Output new fasta
write.table(new.fa
          , "all_sequences_Illumina_Miseq_mod.fasta", row.names=F, quote=F, col.names=F)
