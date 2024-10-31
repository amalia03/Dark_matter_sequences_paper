#!/bin/env Rscript

###kmer analysis using lmjs kmer script (more efficient)
library("scales")
source("kmer_functions.R")


dyn.load("/home/ama/R/R_c_plugins/nuc_kmer_count.so")
args <- commandArgs(trailingOnly=TRUE)


if(length(args)<6){
    stop(paste("You need 6 arguments, the FASTA1 file, the output1 id, FASTA2 file, output2 id, the minimum kmer size (integer) and max kmer size", "\n",
               "Example: ./apply_kmer_any.R fasta_file.fa1 output_id1 fasta_file2 output_id2 5 8") , call.=FALSE)
}

write("Importing files.", stdout())

##Imported data space
##comp.group <- readLines("../meta_symmetry/GEANS/GEANS/120A_head_nless.fa")
#comp.group <- readLines("all_unaligned_300.fa")
#comp.name <- "unaligned"

comp.group <- readLines(args[1])
comp.name <- args[2]
##Last one is an empty string so I am removing that one
#comp2.group <- readLines("all_aligned_300.fa")
#comp2.name <- "aligned"
comp2.group <- readLines(args[3])
comp2.name <- args[4]

comp.seqs <- comp.group[seq(2,length(comp.group), 2)]
comp2.seqs <- comp2.group[seq(2,length(comp2.group), 2)]

#min.k <- args[5]
min.k <- 5
#max.k <- args[6]
max.k <- 8
####kmerification

kmer.l <- function(x){
    lapply(min.k:max.k, function(k.s){
        kmerify(x,k.s,norm=F)
                                        #        kmerify(x,k.s)
    })
}
write(paste0("Start kmerization for kmer size ",min.k,"-",max.k), stdout())

comp.kmr <- kmer.l(comp.seqs)
names(comp.kmr) <- c(min.k:max.k)
comp2.kmr <- kmer.l(comp2.seqs)
names(comp2.kmr) <- c(min.k:max.k)

#kmer.l
#kmerify
comp.kmr <- lapply(comp.kmr, function(x){x[-which(x$counts==max(x$counts)),]})

comp.kmr <- kmer.l(comp.seqs)
names(comp.kmr) <- c(min.k:max.k)
comp2.kmr <- kmer.l(comp2.seqs)
names(comp2.kmr) <- c(min.k:max.k)

#kmer.l
#kmerify
comp.kmr <- lapply(comp.kmr, function(x){x[-which(x$counts==max(x$counts)),]})
comp2.kmr <- lapply(comp2.kmr, function(x){x[-which(x$counts==max(x$counts)),]})


head(comp.kmr[[1]]$counts)
head(comp2.kmr[[1]]$counts)

##For a 4 display
#jpeg(paste0(comp.name,"_vs_", comp2.name,"_for_kmr", min.k,"_to_",max.k,".jpg"), height = 1280, width = 960)
#par(mfrow=c(2,2))
##For a single frame
jpeg(paste0(comp.name,"_vs_", comp2.name,"_for_kmr", min.k,"_to_",max.k,".jpg"), height = 800, width = 800)


par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
sapply(1:length(comp.kmr), function(x){
#    plot(log(comp.kmr[[x]]$counts,10), log(comp2.kmr[[x]]$counts,10), xlab=comp.name, ylab=comp2.name, col=alpha("maroon", .2), pch=19, cex.axis=2, cex.lab=2)
    plot(log(comp.kmr[[x]]$counts,10), log(comp2.kmr[[x]]$counts,10),xlab="",ylab ="", col=alpha("maroon", .2), pch=19, cex.axis=2.5)
    #abline(0,sum(log(comp2.kmr[[x]]$counts,10))/ sum(log(comp.kmr[[x]]$counts,10)),col="darkgrey", lty=2, lwd=2)
    #abline(lm(log(comp.kmr[[x]]$counts,10) ~ log(comp2.kmr[[x]]$counts,10)),col="darkgrey", lty=2, lwd=2)
    ##    lr <- lm(comp.kmr[[x]]$counts~comp2.kmr[[x]]$counts)[[1]][[2]]
#    lr <- summary(lm(comp.kmr[[x]]$counts~comp2.kmr[[x]]$counts))[[9]]
    lr <- round(cor(comp.kmr[[x]]$counts, comp2.kmr[[x]]$counts), 2)
#    p.v <- round(summary(lm(comp.kmr[[x]]$counts~comp2.kmr[[x]]$counts))$coefficients[1,"Pr(>|t|)"],3)
#    p.v <- round(summary(lm(comp.kmr[[x]]$counts~comp2.kmr[[x]]$counts))$coefficients[1,"Pr(>|t|)"],3)
    p.v <- round(as.numeric(cor.test(comp.kmr[[x]]$counts, comp2.kmr[[x]]$counts)[3]), 2)
#    with(par(),text(usr[2]*.95, usr[4]*.95,labels=paste("R2=", round(lr,2)), cex=1.5, font=2))
    legend("topright", c(paste("r=", lr),
                         paste0("p ", ifelse(p.v==0, "<0.001", p.v))),
                         cex=2.4, bty="n", border="white")

})

dev.off()

