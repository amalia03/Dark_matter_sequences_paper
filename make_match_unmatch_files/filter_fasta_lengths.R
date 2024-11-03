#!/bin/R
## This script takes a fasta file as an input and only outputs sequences of a specific length. 

read.fa <- function(f){
    l <- readLines( f )
    id <- l[seq(1, length(l), 2)]
    seq <- l[seq(2, length(l), 2)]
    names(seq) <- id
    seq
}

assembled.fa <- read.fa("all_merged.fa")
read.l <- nchar(assembled.fa)

plot(sort(read.l), type="l")
table.length <- table(read.l)

l <- length(read.l[read.l >= 300])
r.300 <- read.l[read.l >= 300]
plot(sort(r.300), type="l", ylab="Sequence length (bp)")
abline(h=mean(r.300), col="darkgrey", lty=2)
text(100000, 375, "mean")

assembled.300 <- assembled.fa[read.l>=300]
head(cbind(assembled.300, names(assembled.300)))

write.table(assembled.300,"assembled_300.tsv", row.names=T, quote=F, sep="\t", col.names=F)
