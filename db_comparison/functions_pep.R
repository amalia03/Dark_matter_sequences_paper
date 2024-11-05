#!/bin/R
##Functions that are used in the interpro_analysis.R script

##I usually start with this funciton as it is super useful when working on ESS emacs. 
set.display <- function(number){
        local.h <- paste0("localhost:", number,".0")
            Sys.setenv('DISPLAY'=local.h)
        }

# Program to convert decimal number into binary number using recursive function
convert.to.binary <- function(n) {
    if(n > 1) {
        convert.to.binary(as.integer(n/2))
    }
    cat(n %% 2)
}

read.tsv.tables <- function(fn){
    tmp <- lapply( fn, read.table, header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="", fill=TRUE,
                  col.names=c('protid', 'md5', 'len', 'analysis', 'signature', 'description',
                              'start', 'stop', 'score', 'status', 'date', 'accession', 'ipro_des', 'go'))
    do.call(rbind, tmp)
}

## we also want to get the full set of sequences usedin the alignments
read.fa.l <- function(fn){
    tmp <- sapply( fn, function(x){
        lines <- readLines(x)
        lengths <- nchar( lines[ seq(2, length(lines), 2) ])
        names(lengths) <- sub("^>", "", lines[ seq(1, length(lines), 2) ])
        lengths
    })
    tmp2 <- do.call(c, tmp)
        names(tmp2) <- sub("^.+[0-9]+\\.fa\\.([^ ]+).+", "\\1", names(tmp2))
    tmp2
}

##Get top values per query by modular evalue.

top.function <- function(ipro.list, e.thresh=-20){
    lapply(1:length(ipro.list),function(x){
        ipro <- ipro.list[[x]]
        ipro.sc <- tapply(1:nrow(ipro), ipro$protid, function(i){
            i[ order(ipro[i, 'score']) ]
        })
        ipro.sc <- sapply(ipro.sc, function(x){ x[1] })
        ipro.ts <- ipro[ipro.sc,]
        ipro.top <- ipro.ts[ipro.ts$score<10^e.thresh,]
        ipro.top
    })
}

top.function.lite <- function(ipro, e.thresh=-20,ecomp="lower", decrease=F){
        ipro.sc <- tapply(1:nrow(ipro), ipro$protid, function(i){
            i[ order(ipro[i, 'score'], decreasing=decrease) ]
        })
        ipro.sc <- sapply(ipro.sc, function(x){ x[1] })
        ipro.ts <- ipro[ipro.sc,]
        if(ecomp=="lower"){
            ipro.top <- ipro.ts[ipro.ts$score<10^e.thresh,]}
        else if(ecomp=="greater"){
            ipro.top <- ipro.ts[ipro.ts$score>10^e.thresh,]}
        ipro.top
}
