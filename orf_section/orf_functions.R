#Some functions for R. 

read.orfs <- function(f){
    l <- readLines( f )
    id.i <- grep("^>", l)
    id.l <- l[id.i]
    sw <- strsplit( id.l, '\t' )
    coords <- sapply(sw, function(x){ as.numeric(x[2:4]) })
    ids <- sapply(sw, function(x){ x[1] })
    colnames(coords) <- ids
    t(coords)
}

##  Assumes unfolded fasta
read.fa <- function(f){
    l <- readLines( f )
    id <- l[seq(1, length(l), 2)]
    seq <- l[seq(2, length(l), 2)]
    names(seq) <- id
    seq
}
