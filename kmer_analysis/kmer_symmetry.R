#!/bin/env Rscript
#Make sure you have the kmer_functions.R on the same directory or rewrite the path. 
source("kmer_functions.R")
#Make sure you have R_c_plugins (link in the README file) in your R directory.

dyn.load("../R_c_plugins/nuc_kmer_count.so")

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("You need 4 arguments, fasta file, output name, col1 and col2 \n", call.=FALSE)
}


fasta.file <- readLines(args[1])
fa.seqs <- fasta.file[seq(2,length(fasta.file), 2)]

jpeg(paste0(args[2], "_symmetry.jpg"),  width = 1024, height = 768)
par(mfrow=c(2,3))
sapply(4:9, function(k.s){
        write(paste("Starting analysis with kmer size =", k.s, "at", Sys.time()), stdout())
            symmetrify.dens(fa.seqs,k.s, maint=paste0("kmer and reverse complement at size= ",k.s, args[2]), col1=args[3], col2=args[4], pix
=100)
        })
dev.off()
