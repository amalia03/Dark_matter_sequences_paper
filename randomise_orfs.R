####An r script that creates ten randomized replicates based on the length of each ORF sequence, using the weights for each file. 

if(lengths(args)==0){
    stop("You need 3 entries, the nucleotide frequency file, a file containing the sequence lengths and output name. ")
}

nucl<- read.table(args[1], sep="\t", header=F, stringsAsFactors=F)
###         nucl<- read.table("300bp_nucl.tsv", sep="\t", header=F, stringsAsFactors=F)
nucl

colnames(nucl) <- c("percent_A","percent_C","percent_T","percent_G")
)

len<- read.table(args[2], sep="\t", header=F, stringsAsFactors=F)
colnames(len)<- c("id", "seq_len")
                                            
write(paste("Starting randomisation process at", Sys.time()), stdout())

random.seqs <- data.frame(id=c(paste0("id_",1:(nrow(len)*10))),
                          seq=c(sapply(
                              1:(nrow(len)),
                              function(i){
                                  replicate(10,paste(sample(c("A","C","G","T"), len$seq_len[i], replace=TRUE, prob = c(nucl$percent_A, nucl$percent_C, nucl$percent_G, nucl$percent_T)),collapse=""))
                              }
                          )
                          )
                          )

##
random.sample <- random.seqs[sample(nrow(random.seqs), nrow(random.seqs)/5),]

dim(random.seqs)
dim(random.sample)

write.table(random.sample, paste0(args[3], "_random_lengths.tsv"),
            col.names=FALSE, row.names=F, quote=F, sep="\t")
