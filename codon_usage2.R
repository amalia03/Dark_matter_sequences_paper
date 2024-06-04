## consider codon usage from different subsets of the data:
                                                                                                     
g.code <- read.table("./genetic_code.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)

codon.usage.f <- list.files(".", pattern="codon.tsv", full.names=TRUE)
names(codon.usage.f) <- c('al', 'random', 'zebrafish', 'unaligned')
codon.usage.labels <- c('Aligned', "Random", "Zebrafish genome", "Unaligned")
names(codon.usage.labels) <- names(codon.usage.f)
codon.usage.f


codon.usage <- lapply( codon.usage.f, read.table, sep="\t", header=TRUE, stringsAsFactors=FALSE )

## we also want to incorporate the blast data.
## to select subsets of the sequences
rand.bl <- read.table( "../assembled_300.bl", stringsAsFactors=FALSE, sep="\t",
                      quote="", comment.char="",
                      col.names=c('query', 'subject', 'stitle', 'qlen', 'slen',
                                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'score', 'length', 'pident', 'nident','gapopen', 'gaps', 'qcovs'))

seq.subsets.str <- c(mito="mitoch", rrna="rrna", human="sapiens", nanno="nannochlor")

seq.subsets <- lapply(seq.subsets.str, function(x){
    codon.usage$al$id %in% rand.bl[ grep(x, rand.bl$stitle, ignore.case=TRUE), 'query']
})

## let us make filtered version of the data for the aligned sequences
b <- with( seq.subsets, !(mito | rrna | human) )
head(codon.usage$aligned.f)
codon.usage <- c( 'aligned.f'=list(codon.usage$al[b, ]), codon.usage )
codon.usage.labels <- c('aligned.f'='Aligned filtered', codon.usage.labels)
names(codon.usage) <- c("aligned.f", "al", "random", "unaligned","zebrafish")
## First lets have a quick look at the distribution of ORF coverages.
breaks <- seq(0, 1, 0.01)
min.l <- 300

pdf("coding_potential_400_ama.pdf", width=12, height=8)
par(mfrow=c(2, 3))
invisible( lapply(names(codon.usage), function(nm){
    cod <- codon.usage[[nm]]
    b1 <- with(cod, nl > min.l & rank == 1 & frame < 0)
    b2 <- with(cod, nl > min.l & rank == 1 & frame > 0)
    h1 <- with( cod, hist( cov[ b1 ], breaks=breaks, plot=FALSE ))
    h2 <- with( cod, hist( cov[ b2 ], breaks=breaks, plot=FALSE ))
    plot(1,1, type='n', xlab='Coverage', ylab='Count', main=codon.usage.labels[nm],
         xlim=c(0,1), ylim=range(c(h1$counts, h2$counts)), cex.lab=1.5)
    l <- length(breaks)
    with(cod, legend('topleft',
                     legend=c(sprintf("%.1f %%", 100 * sum(cov[b1] == 1) / sum(b1)),
                              sprintf("%.1f %%", 100 * sum(cov[b2] == 1) / sum(b2))),
                     text.col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))))
    with(h1, rect( breaks[-l], 0, breaks[-1], counts, col=rgb(1, 0, 0, 0.5), border=NA))
    with(h2, rect( breaks[-l], 0, breaks[-1], counts, col=rgb(0, 0, 1, 0.5), border=NA))
}))
dev.off()

## plot as percentage..
pdf("coding_potential_400_ama.pdf", width=12, height=8)
breaks <- seq(0, 1, 0.01)
min.l <- 300
par(mfrow=c(2, 3))
invisible( lapply(names(codon.usage), function(nm){
    cod <- codon.usage[[nm]]
    b1 <- with(cod, nl > min.l & rank == 1 & frame < 0)
    b2 <- with(cod, nl > min.l & rank == 1 & frame > 0)
    h1 <- with( cod, hist( cov[ b1 ], breaks=breaks, plot=FALSE ))
    h2 <- with( cod, hist( cov[ b2 ], breaks=breaks, plot=FALSE ))
    plot(1,1, type='n', xlab='Coverage', ylab='Density', main=codon.usage.labels[nm],
         xlim=c(0,1), ylim=range(c(h1$counts/sum(h1$counts), h2$counts/sum(h2$counts))), cex.lab=1.2, cex.axis=1.3)
    l <- length(breaks)
    with(cod, legend('topleft',
                     legend=c(sprintf("-ve ORF: %.1f %%", 100 * sum(cov[b1] == 1) / sum(b1)),
                              sprintf("+ve ORF: %.1f %%", 100 * sum(cov[b2] == 1) / sum(b2))),
                     text.col=c(rgb(1,0,0,0.7), rgb(0,0,1,0.7)), cex=1.1, font=4, bg="white"))
    with(h1, rect( breaks[-l], 0, breaks[-1], counts/sum(counts), col=rgb(1, 0, 0, 0.5), border=NA))
    with(h2, rect( breaks[-l], 0, breaks[-1], counts/sum(counts), col=rgb(0, 0, 1, 0.5), border=NA))
}))
dev.off()

## do also, but with aligned reads filtered so that they do not include
## mitochondrial or other contaminants defined above.
##b <- with( seq.subsets, !(mito | rrna | human | nanno) )


## This now at least makes sense. For sequences longer than 400
## aligned:    48.5%
## unaligned:  32.9%
## have complete Open reading frames
## Actually quite a big difference, but:
## Suggests that at least 60% of unaligned sequences are derived
## from the appropriate RNA strand.

## take a table with the appropriate counts and the codons to be counted
get.codon.freq <- function(df, codons){
    counts <- colSums( df[,codons] )
}

## codon.usage frequencies for all reads above min.l
## Corresponds to :
## aligned.f   aligned        hg      rand      unal 
##     30144     73050     36666    468780    154524 
## (minimum length requirement used since it makes it
## more likely that we pick the real ORF
codon.usage.f1 <- lapply( codon.usage, function(x){
    b <- x$nl >= min.l
    tmp <- mapply( function(is.fwd, rank){
        get.codon.freq( x[x$frame * is.fwd > 0 & x$rank == rank & b, ], rownames(g.code) )},
        c(1,1,1,-1,-1,-1), c(1:3,1:3))
    colnames(tmp) <- paste(c('f', 'f', 'f', 'r', 'r', 'r'), c(1:3,1:3), sep="_")
    tmp
})

## codon.usage frequencies for long complete coverage sequences
codon.usage.f2 <- lapply( codon.usage, function(x){
    b <- x$nl >= min.l
    cov.b <- list(x$cov == 1 & b, x$cov > 0.5 & b, x$cov > 0.5 & b)
    tmp <- mapply( function(is.fwd, rank, b){
        get.codon.freq( x[x$frame * is.fwd > 0 & x$rank == rank & b, ], rownames(g.code) )},
        c(1,1,1,-1,-1,-1), c(1:3,1:3), c(cov.b, cov.b))
    colnames(tmp) <- paste(c('f', 'f', 'f', 'r', 'r', 'r'), c(1:3,1:3), sep="_")
    tmp
})

## codon usage table (cu)
## genetic code (g.code)
codon.bias <- function(cu, code){
    cu <- t(cu)
    cb <- tapply(1:ncol(cu), g.code[ colnames(cu), 'aa' ], function(i){
        rs <- rowSums( cu[,i, drop=FALSE] )
        rs <- ifelse(rs == 0, 1, rs)
        t(cu[,i, drop=FALSE] / rs)
    })
}
        

## Also calculate codon bias
codon.usage.bs.1 <- lapply(codon.usage.f1, codon.bias, code=g.code )
codon.usage.bs.2 <- lapply(codon.usage.f2, codon.bias, code=g.code )

## we are interested in the codon usage correlation primarily between unaligned
## and aligned sequences:

cu.cor.1.1 <- with( codon.usage.f1, cor(cbind(aligned, unal)) )
cu.cor.1.2 <- with( codon.usage.f1, cor(cbind(aligned.f, unal)) )
cu.cor.1.3 <- with( codon.usage.f1, cor(cbind(aligned, aligned.f)) )

cu.cor.2.1 <- with( codon.usage.f2, cor(cbind(aligned, unal)) )
cu.cor.2.2 <- with( codon.usage.f2, cor(cbind(aligned.f, unal)) )
cu.cor.2.3 <- with( codon.usage.f2, cor(cbind(aligned, aligned.f)) )

cb.cor.1.1 <- with(codon.usage.bs.1, mapply(function(x, y){ list(cor(cbind(x,y))) }, aligned, unal))

cb.cor.2.1 <- with(codon.usage.bs.2, mapply(function(x, y){ list(cor(cbind(x,y))) }, aligned, unal))
cb.cor.2.2 <- with(codon.usage.bs.2, mapply(function(x, y){ list(cor(cbind(x,y))) }, aligned.f, unal))
cb.cor.2.3 <- with(codon.usage.bs.2, mapply(function(x, y){ list(cor(cbind(x,y))) }, aligned, aligned.f))
## these are very confusing and I don't really know how to make use of them.


plot.cor <- function(x, cex=1.25){
    image(x=1:ncol(x), y=1:nrow(x), x, xaxt='n', yaxt='n', xlab='', ylab='')
    mtext( colnames(x), 1, at=1:ncol(x), cex=cex, line=1 )
    mtext( rownames(x), 2, at=1:nrow(x), cex=cex, line=1 )
}

plot.cor.rows <- function(x, rows, col=1:length(rows), xlab='ORF', ylab='Correlation',
                          lwd=1){
    plot(1:ncol(x), x[rows[1], ], type='n', xlab=xlab, ylab=ylab, xaxt='n')
    for(i in 1:length(rows)){
        lines(1:ncol(x), x[rows[i],], col=col[1 + (i-1) %% length(col)], lwd=lwd)
    }
    axis(side=1, at=1:ncol(x), labels=colnames(x))
}

par(mfrow=c(2,3))
plot.cor(cu.cor.1.1, cex=1)
plot.cor(cu.cor.1.2, cex=1)
plot.cor(cu.cor.1.3, cex=1)

plot.cor.rows( cu.cor.1.1, c('r_1'), lwd=3 )
plot.cor.rows( cu.cor.1.2, c('r_1'), lwd=3 )
plot.cor.rows( cu.cor.1.3, c('r_1'), lwd=3 )

par(mfrow=c(2,3))
plot.cor(cu.cor.2.1, cex=1)
plot.cor(cu.cor.2.2, cex=1)
plot.cor(cu.cor.2.3, cex=1)

plot.cor.rows( cu.cor.2.1, c('r_1'), lwd=3 )
plot.cor.rows( cu.cor.2.2, c('r_1'), lwd=3 )
plot.cor.rows( cu.cor.2.3, c('r_1'), lwd=3 )

### These plots suggest that there is a difference in the codon usage of those sequences
### that had alignments to nr/nt, and that this difference is not affected much by the
### removal of mitochondrial sequences

### Let us see how robust this difference is by bootstrapping the codon
### usage frequencies:

## cu is a codon.usage table
## orientation = -1 or +1
## rank 1,2, or 3
## bs.size is a proportion of the rows that pass the initial filter.
bootstrap.frequency <- function(cu, orientation, rank, min.nl, min.cov, bs.size, bs.n){
    b <- cu$frame * orientation > 0 & cu$nl >= min.l & cu$cov >= min.cov & cu$rank == rank
    b.i <- which(b)
    bs.sn <- as.integer( length(b.i) * bs.size ) ## the sample number used
    total.f <- get.codon.freq(cu[b, ], rownames(g.code))
    bs.f <- sapply( 1:bs.n, function(i){
        get.codon.freq(cu[ sample(b.i, bs.sn, replace=TRUE), ], rownames(g.code))
    })
    codon.cor <- cor( total.f, bs.f )
    list('i'=b.i, 'full'=total.f, 'bs'=bs.f, 'cor'=codon.cor)
}

codon.usage.f.bs.1 <- lapply( codon.usage, function(x){
    ## only look at r_1 (reverse best ORF)
    ## minimum length of 300, full coverage
    bootstrap.frequency( x, -1, 1, 300, 1, 0.5, 100 )
})

codon.usage.f.bs.2 <- lapply( codon.usage, function(x){
    ## only look at r_1 (reverse best ORF)
    ## minimum length of 300, coverage >= 0.75
    bootstrap.frequency( x, -1, 1, 300, 0.75, 0.5, 100 )
})

par(mfrow=c(2,3))
for(nm in names(codon.usage.f.bs.1)){
    hist(codon.usage.f.bs.1[[nm]]$cor, main=nm)
}

par(mfrow=c(2,3))
for(nm in names(codon.usage.f.bs.2)){
    hist(codon.usage.f.bs.2[[nm]]$cor, main=nm)
}

## internally we have (as one would expect), very, very high correlations
## Interestingly, these are highest for the unaligned sequences.
## Let us look at the cross correlation of aligned.f and unal
## These will presumably be very similar:

al.unal.bs.cor.1 <- with( codon.usage.f.bs.1, cor( aligned$bs, unal$bs ) )
al.unal.bs.cor.2 <- with( codon.usage.f.bs.2, cor( aligned$bs, unal$bs ) )

al.f.unal.bs.cor.1 <- with( codon.usage.f.bs.1, cor( aligned.f$bs, unal$bs ) )
al.f.unal.bs.cor.2 <- with( codon.usage.f.bs.2, cor( aligned.f$bs, unal$bs ) )

hist(al.unal.bs.cor.1)
hist(al.unal.bs.cor.2)
hist(al.f.unal.bs.cor.1)
hist(al.f.unal.bs.cor.2)
## these are very similar: all complete ORFs gives 0.76, and with at least 75%
## coverage we have 0.78. 

## So there is something that is different about the unaligned ORFs as a
## set. Is there a subset of the aligned ORFs that looks more similar
## to the unaligned ORFs
##
## Read in taxonomy information for the blast data and try to look at
## distinct subsets.

rand.bl.tax <- read.table("../../../blast_taxonomy/random_miseq_seqs_vs_nt2.taxonomy",
                          stringsAsFactors=FALSE, sep="\t", quote="", comment.char="")
colnames(rand.bl.tax) <- c('subject', 'tax.id', 'taxon', 'nodes.id', 'nodes.name', 'query', 'evalue', 'score', 'length', 'pident',
                           'nident')

## Somehow:
length(unique(rand.bl.tax$query))
## [1] 379449
## is larger than:
length(unique(codon.usage$aligned$id))
## [1] 358046
## I need to work out why that is.
head( setdiff( rand.bl.tax$query, unique( codon.usage$aligned$id )) )
## [1] "M02443:118:000000000-AMUKU:1:1101:19240:1804"
## [2] "M02443:118:000000000-AMUKU:1:1101:8143:1984" 
## [3] "M02443:118:000000000-AMUKU:1:1101:11836:3333"
## [4] "M02443:118:000000000-AMUKU:1:1101:9814:3682" 
## [5] "M02443:118:000000000-AMUKU:1:1101:10636:3837"
## [6] "M02443:118:000000000-AMUKU:1:1101:24732:3868"
##
## Looks like it may be related to these sequences being too short. At least the
## first one of the R1s has a length of 41 after trimming. That probably means
## that the mergins software did not merge it.


tax.terms <- unlist( strsplit( rand.bl.tax$nodes.name, "," ))
length(tax.terms)
## [1] 6526474
tax.terms.t <- sort( table(tax.terms), decreasing=TRUE )

## calculate codon usage frequencies for the top set of terms for
## the aligned data set. Only for the best reverse ORF longer than
## 300 bp and with complete coverage:

min.l <- 300;
## I should probably return some more information, like the indices of
## the matching sequences. But lets first see if it works.
al.codon.usage.1 <- lapply( (names(tax.terms.t)[1:200]), function(term){
    cu <- codon.usage$aligned.f
    b <- cu$nl >= min.l & cu$cov == 1 & cu$rank == 1 & cu$frame < 0
    qid <- rand.bl.tax$query[ grep(term, rand.bl.tax$nodes.name) ]
    b <- b & cu$id %in% qid
    get.codon.freq( cu[b, ], rownames(g.code))
})
names(al.codon.usage.1) <- (names(tax.terms.t)[1:200])

al.codon.usage.1.unal.cor <- cor( codon.usage.f2$unal[,'r_1'], do.call(cbind, al.codon.usage.1))
barplot( al.codon.usage.1.unal.cor)

## al.codon.usage.1.unal.cor is a matrix with one row.. hence sort removes the column names
## but we can:
al.codon.usage.1.unal.cor[, order(al.codon.usage.1.unal.cor)]
## to see the patter.

## The best correlation is for Alveolata
i1 <- grep("Alveolata", rand.bl.tax$nodes.name)
## 6515 matches

head( rand.bl.tax[i1, c('subject', 'taxon', 'evalue', 'score', 'length', 'pident') ], n=10 )

## which is nice. However, the second best one is for Apocrita which are
## insects. Which is not nice (as we shouldn't see them here)
i2 <- grep("Apocrita", rand.bl.tax$nodes.name)
## 2504 matches
head( rand.bl.tax[i2, c('subject', 'taxon', 'evalue', 'score', 'length', 'pident') ], n=10 )

## compare the qualities of the scores:
with( rand.bl.tax, plot( length[i1], pident[i1] ))
with( rand.bl.tax, plot( length[i2], pident[i2] ))

with( rand.bl.tax, plot( pident[i1], -log10(evalue[i1]) ))
with( rand.bl.tax, plot( pident[i2], -log10(evalue[i2]) ))

## these suggest that maybe we can simply set a threshold of 95% and see how that affects
## the correlations and numbers that we get.

min.l <- 300;
min.pid <- 95
## I should probably return some more information, like the indices of
## the matching sequences. But lets first see if it works.
al.codon.usage.2 <- lapply( (names(tax.terms.t)[1:200]), function(term){
    cu <- codon.usage$aligned.f
    b <- cu$nl >= min.l & cu$cov == 1 & cu$rank == 1 & cu$frame < 0
    bl.b <- rand.bl.tax$pident >= min.pid
    qid <- (rand.bl.tax$query[bl.b])[ grep(term, rand.bl.tax$nodes.name[bl.b]) ]
    b <- b & cu$id %in% qid
    get.codon.freq( cu[b, ], rownames(g.code))
})
names(al.codon.usage.2) <- (names(tax.terms.t)[1:200])

al.codon.usage.2.unal.cor <- cor( codon.usage.f2$unal[,'r_1'], do.call(cbind, al.codon.usage.2))
barplot( al.codon.usage.2.unal.cor)

al.codon.usage.2.unal.cor[, order(al.codon.usage.2.unal.cor)]
## a lot of NAs due to no levels, but, the best are Foraminifera, Retaria, and Rhizaria
## which are all microplankton of sorts (well these three groups seem to be completely
## overlapping (with foraminerea being the most distinct group)
## but this is then followed by:
## Ovalentaria (teleost), viruses, eudicotyledons (plants)
## But for this the codon counts are really too small. 

min.l <- 300;
min.pid <- 90
max.evalue <- 1e-20
## I should probably return some more information, like the indices of
## the matching sequences. But lets first see if it works.
al.codon.usage.3 <- lapply( (names(tax.terms.t)[1:200]), function(term){
    cu <- codon.usage$aligned.f
    b <- cu$nl >= min.l & cu$cov == 1 & cu$rank == 1 & cu$frame < 0
    bl.b <- rand.bl.tax$pident >= min.pid & rand.bl.tax$evalue >= max.evalue
    qid <- (rand.bl.tax$query[bl.b])[ grep(term, rand.bl.tax$nodes.name[bl.b]) ]
    b <- b & cu$id %in% qid
    get.codon.freq( cu[b, ], rownames(g.code))
})
names(al.codon.usage.3) <- (names(tax.terms.t)[1:200])

al.codon.usage.3.unal.cor <- cor( codon.usage.f2$unal[,'r_1'], do.call(cbind, al.codon.usage.3))
barplot( al.codon.usage.3.unal.cor)

al.codon.usage.3.unal.cor[, order(al.codon.usage.3.unal.cor)]
## this on the other hand, gives mainly a load of teleosts... One wonders why.

## In a perverse way it makes sense that we get a load of things that do not make sense
## Because the matches to species that do not make sense are more likely to have arisen
## from contaminations. That is, we suspect that matches to rice would arise from contamination
## of rice with somethink like water mould;
## But seeing teleosts is not so good, since we know that a large part of the sequences were
## derived from the lumpsucker; of course these may be included in both of the data sets here.

## I have now obtained codon counts for the subject (i.e. nr/nt) sequences. Many of the
## subject ids were actually genome sequences containing many ORFs. I have obtained counts
## based on the full set of ORFs for these:

nt.codon.u.f <- list.files("../../../mapped_gb/mp_dir", pattern="cu_2_.+tsv", full.names=TRUE)

nt.codon.u <- lapply( nt.codon.u.f, read.table, sep="\t", header=TRUE, stringsAsFactors=FALSE )
nt.codon.u <- do.call(rbind, nt.codon.u)

dim(nt.codon.u)
## [1] 62178    71
sum( nt.codon.u$accession %in% rand.bl.tax$subject )
## 62019
## so what's missing?

with(nt.codon.u, head( accession[ ! accession %in% rand.bl.tax$subject ] ))
## [1] "PFATUBA"    "DDU02283"   "LEIHSP70G"  "MLGMOCUMAA" "EIMDEVGEND"
## [6] "CBRG03N05" 
##
## I unfortunately used the LOCUS line to obtain the accession id to make things
## simpler (since the accession may contain several different versions)
## Unfortunately it turns out that the id on the LOCUS lines sometimes differs
## from the accession. And that is probably the case for these. I do not wish to
## rerun the extraction, so for now, I will live with it by simply removing those
## that do not work.

nt.codon.u <- nt.codon.u[ nt.codon.u$accession %in% rand.bl.tax$subject, ]
dim(nt.codon.u)
## [1] 62019    71
##
nt.codon.u.tax <- with(rand.bl.tax, nodes.name[ match( nt.codon.u$accession, subject ) ])

## we can now try to calculate codon usage for different classes:
nt.codon.usage.1 <- t(sapply( (names(tax.terms.t)[1:200]), function(term){
    b <- grepl( term, nt.codon.u.tax )
    get.codon.freq( nt.codon.u[b,], rownames(g.code) )
}))
rownames(nt.codon.usage.1) <- (names(tax.terms.t)[1:200])

## do these have any obvious patterns?
## note that this will include mitochondrial and other sequences to confuse issues, but
## we can nevertheless consider if there is any particular pattern.

nt.codon.pca.1 <- prcomp(nt.codon.usage.1, scale=TRUE)

plot(nt.codon.pca.1)  ## almost everything is in the first dimension.
with(nt.codon.pca.1, plot(x[,1], x[,2], cex=2))
with(nt.codon.pca.1, identify(x[,1], x[,2], labels=rownames(x)) )
## The division here is taken up almost entirely by bacterial groups
## lets redo the pca excluding bacterial labels:

b <- !(grepl("bacteria", rownames(nt.codon.usage.1), ignore.case=TRUE ))
sum(b) ## 195
nt.codon.pca.1.1 <- prcomp( nt.codon.usage.1[b,], scale=TRUE )

plot(nt.codon.pca.1.1)  ## still almost everything is in the first dimension.
with(nt.codon.pca.1.1, plot(x[,1], x[,2], cex=2))
with(nt.codon.pca.1.1, identify(x[,1], x[,2], labels=rownames(x)) )

## zoom in a bit on the plot
with(nt.codon.pca.1.1, plot(x[,1], x[,2], cex=2, xlim=c(-20,10), ylim=c(-7,10)))
with(nt.codon.pca.1.1, identify(x[,1], x[,2], labels=rownames(x)) )
with(nt.codon.pca.1.1, points( x['Foraminifera',1], x['Foraminifera', 2], col='red', pch=19 ))

## and much, much more to see the bunched up stuff:
with(nt.codon.pca.1.1, plot(x[,1], x[,2], cex=2, xlim=c(-3,5), ylim=c(-3,3)))
with(nt.codon.pca.1.1, identify(x[,1], x[,2], labels=rownames(x)) )
with(nt.codon.pca.1.1, points( x['Foraminifera',1], x['Foraminifera', 2], col='red', pch=19 ))

## and again, with a much, much smaller space:
with(nt.codon.pca.1.1, plot(x[,1], x[,2], cex=2, xlim=c(2.75,3.5), ylim=c(-0.5,0)))
with(nt.codon.pca.1.1, identify(x[,1], x[,2], labels=rownames(x)) )
with(nt.codon.pca.1.1, points( x['Foraminifera',1], x['Foraminifera', 2], col='red', pch=19 ))

## Do the same, but this time try to avoid including mitochondrial sequences

## we can now try to calculate codon usage for different classes:
rand.bl.accession <- sub("[^|]*\\|([^|]+)\\|?.*", "\\1", rand.bl$subject )
rand.bl.accession <- sub("\\.[0-9]+", "", rand.bl.accession)
sum( nt.codon.u$accession %in% rand.bl.accession )
## 62019 .. finally all of them match. 

i <- match(nt.codon.u$accession, rand.bl.accession )
mito.b <- grepl("mitoc", rand.bl$stitle[i], ignore.case=TRUE) ## 3402. not so many
nt.codon.usage.2 <- t(sapply( (names(tax.terms.t)[1:200]), function(term){
    b <- (!mito.b) & grepl( term, nt.codon.u.tax )
    get.codon.freq( nt.codon.u[b,], rownames(g.code) )
}))
rownames(nt.codon.usage.2) <- (names(tax.terms.t)[1:200])

## and let us also do this but with only entries for single sequences
## to avoid the large number of genomic fragments:
i <- match(nt.codon.u$accession, rand.bl.accession )
mito.b <- grepl("mitoc", rand.bl$stitle[i], ignore.case=TRUE) ## 3402. not so many
nt.codon.usage.3 <- t(sapply( (names(tax.terms.t)[1:200]), function(term){
    b <- (!mito.b) & nt.codon.u$nORF == 1 & grepl( term, nt.codon.u.tax )
    get.codon.freq( nt.codon.u[b,], rownames(g.code) )
}))
rownames(nt.codon.usage.3) <- (names(tax.terms.t)[1:200])



b <- !(grepl("bacteria", rownames(nt.codon.usage.2), ignore.case=TRUE ))
sum(b) ## 195
nt.codon.pca.2.1 <- prcomp( nt.codon.usage.2[b,], scale=TRUE )

plot(nt.codon.pca.2.1)  ## still almost everything is in the first dimension.
with(nt.codon.pca.2.1, plot(x[,1], x[,2], cex=2))
with(nt.codon.pca.2.1, identify(x[,1], x[,2], labels=rownames(x)) )

with(nt.codon.pca.2.1,  plot(x[,1], x[,2], xlim=c(-1.5,3.75), ylim=c(-2,2), cex=2))
with(nt.codon.pca.2.1, identify(x[,1], x[,2], labels=rownames(x)) )

with(nt.codon.pca.2.1,  plot(x[,1], x[,2], xlim=c(2.6,3.3), ylim=c(-0.4,-0.1), cex=2))
with(nt.codon.pca.2.1, identify(x[,1], x[,2], labels=rownames(x)) )

## this also looks rather messy, but with some not completely unreasonable patterns.
## first compare these to codon.usage.f2 (long) and codon.usage.f2 (long and complete)

nt.codon.cor.1.1 <- cor( codon.usage.f1$unal[,'r_1'], t(nt.codon.usage.1) )
## and that gives us viruses and arthropods as the best scoring things.
nt.codon.cor.1.2 <- cor( codon.usage.f1$unal[,'r_1'], t(nt.codon.usage.2) )
## a matrix, so we must
nt.codon.cor.1.2[ ,order(nt.codon.cor.1.2) ]

## let us also consider:
nt.codon.cor.2.1 <- cor( codon.usage.f2$unal[,'r_1'], t(nt.codon.usage.1) )
## and that gives us viruses and arthropods as the best scoring things.
nt.codon.cor.2.2 <- cor( codon.usage.f2$unal[,'r_1'], t(nt.codon.usage.2) )

nt.codon.cor.2.2[ ,order(nt.codon.cor.2.2) ]

## let us also consider:
nt.codon.cor.1.3 <- cor( codon.usage.f1$unal[,'r_1'], t(nt.codon.usage.3) )
## and that gives us viruses and arthropods as the best scoring things.
nt.codon.cor.2.3 <- cor( codon.usage.f2$unal[,'r_1'], t(nt.codon.usage.3) )

nt.codon.cor.1.3[ ,order(nt.codon.cor.1.3) ]
nt.codon.cor.2.3[ ,order(nt.codon.cor.2.3) ]

### it's all a mess. Let us consider the correlation between the query and
### subject codon correlations:

qs.cor.1.1 <- cor( do.call(cbind, al.codon.usage.1), t(nt.codon.usage.1) )
qs.cor.1.3 <- cor( do.call(cbind, al.codon.usage.1), t(nt.codon.usage.3) )
qs.cor.2.3 <- cor( do.call(cbind, al.codon.usage.2), t(nt.codon.usage.3) )

o <- order(sapply(1:nrow(qs.cor.1.3), function(i){ max( qs.cor.1.3[i,-i], na.rm=TRUE)}), decreasing=TRUE)

qs.cor.1.1[o[1:10], o[1:10]]
qs.cor.1.3[o[1:10], o[1:10]]

qs.cor.1.1[o[1:10], o[1:10]]
qs.cor.1.3[o[1:10], o[1:10]]

qs.cor.2.3[o[1:10], o[1:10]]
## Terrible??

## this really doesn't make that much sense. Give up on it for now.
## but let us look at:

plot( al.codon.usage.1[['Alveolata']], nt.codon.usage.1['Alveolata',] )
plot( al.codon.usage.1[['Alveolata']], nt.codon.usage.2['Alveolata',] )

## The third one here make much more sense than the others.
plot( al.codon.usage.1[['Alveolata']], nt.codon.usage.3['Alveolata',] )
plot( al.codon.usage.1[['Alveolata']], nt.codon.usage.3['Apocrita',] )
plot( al.codon.usage.1[['Alveolata']], nt.codon.usage.3['Foraminifera',] )

plot( al.codon.usage.1[['Foraminifera']], nt.codon.usage.3['Foraminifera',] )
plot( al.codon.usage.2[['Foraminifera']], nt.codon.usage.3['Foraminifera',] )
## the latter does not look good, but that maybe simply because of the smaller
## numbers of counts for al.codon.usage.2 (max 60, compared to more then 2500)

plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Chelicerata',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Arachnida',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Apocrita',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Hymenoptera',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Aconoidasida',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Cnidaria',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Nematoda',] )

plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Haemosporida',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Plasmodiidae',] )
plot( codon.usage.f2$unal[,'r_1'], nt.codon.usage.3['Plasmodium',] )

plot( codon.usage.f1$unal[,'r_1'], nt.codon.usage.3['Aconoidasida',] )
plot( codon.usage.f1$unal[,'r_1'], nt.codon.usage.3['Chelicerata',] )
plot( codon.usage.f1$unal[,'r_1'], nt.codon.usage.3['Haemosporida',] )
plot( codon.usage.f1$unal[,'r_1'], nt.codon.usage.3['Plasmodiidae',] )
identify( codon.usage.f1$unal[,'r_1'], nt.codon.usage.3['Plasmodiidae',], colnames(nt.codon.usage.3) )


with(codon.usage.f2, plot(unal[,'r_1'] / sum(unal[,'r_1']), type='l'))
lines(nt.codon.usage.3['Chelicerata',] / sum(nt.codon.usage.3['Chelicerata',]), col='red')
lines(nt.codon.usage.3['Arachnida',] / sum(nt.codon.usage.3['Arachnida',]), col='blue')
lines(nt.codon.usage.3['Aconoidasida',] / sum(nt.codon.usage.3['Aconoidasida',]), col='green')
lines(nt.codon.usage.3['Nematoda',] / sum(nt.codon.usage.3['Nematoda',]), col='brown', lwd=2)
lines(nt.codon.usage.3['Cnidaria',] / sum(nt.codon.usage.3['Cnidaria',]), col='cyan', lwd=2)
lines(nt.codon.usage.3['Lophotrochozoa',] / sum(nt.codon.usage.3['Lophotrochozoa',]), col='gold', lwd=2, lty=2)

## cu is codon usage
## taxons a list of taxons (giving the rows)
## gcode = g.code 3 column data.frame with
## the codon as the rowname
plot.by.aa <- function(cu, taxons, gcode, nt.height=0.02, lwd=1){
    cu <- cu[taxons, ]
    cu <- cu / rowSums( cu )
    ## define an order that we wish to use
    o <- unlist(tapply( 1:nrow(gcode), gcode[,'aa'], eval ));
    ## in order to indicate the codons on the plot
    ## makes a matrix with 3 rows and 64 columns
    codons <- sapply( strsplit( rownames(gcode)[o], '' ), eval )
    aa <- gcode[o, 'aa']
    aa.l <- 0.5 + c(0, which(aa[-1] != aa[-length(aa)]), length(aa))
    min.y <- min(cu) - (3 * nt.height * diff(range(cu)))
    plot( 1, xlim=c(1, length(o)), ylim=c(min.y, max(cu)), xlab='Codon', ylab='Frequency', xaxt='n', xaxs='i')
    usr <- par('usr')
    h <- usr[4]-usr[3]
    x <- 1:length(o)
    rect( aa.l[-length(aa.l)], usr[3], aa.l[-1], usr[4], col=c(rgb(0.8,0.8,0.8), rgb(1,1,1)), border=NA)
    text( x, usr[4], aa, pos=1 )
    ## draw the nucleotide boxes, just make the black for now
    nt.cols <- c('A'='red', 'C'='green', 'G'='blue', 'T'='cyan')
    for(i in 1:nrow(codons)){
        y2 <- usr[3] + i * (h*nt.height)
        y1 <- y2 - (h*nt.height)
        rect( x-0.5, y1, x+0.5, y2, col=nt.cols[ codons[1 + nrow(codons) - i,] ] )
    }
    for(i in 1:nrow(cu))
        lines( x, cu[i,o], col=i, lwd=lwd, type='b' )
    legend('topleft', taxons, text.col=1:length(taxons), cex=2)
}

## 1: all codons
## 2: non-mitochondrial codons
## 3: non-mitochondrial codons from loci with a single CDS
tmp.tax <- c('unal', 'Chelicerata', 'Arachnida', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa')
tmp.cu <- list()
tmp.cu[[1]] <- rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.1)
tmp.cu[[2]] <- rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.2)
tmp.cu[[3]] <- rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.3)

tmp.tax <- c('unal', 'Chelicerata', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa', 'Aves')
tmp.tax <- c('unal', 'Chelicerata', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa', 'Amoebozoa')
tmp.tax <- c('unal', 'Chelicerata', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa', 'Polychaeta')
par(mfrow=c(3,1))
lwd=3
plot.by.aa( tmp.cu[[1]], tmp.tax, g.code, lwd=lwd )
plot.by.aa( tmp.cu[[2]], tmp.tax, g.code, lwd=lwd )
plot.by.aa( tmp.cu[[3]], tmp.tax, g.code, lwd=lwd )

plot.by.aa( rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.1), c('unal', 'Chelicerata', 'Arachnida', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa'), g.code, lwd=lwd )
plot.by.aa( rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.2), c('unal', 'Chelicerata', 'Arachnida', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa'), g.code, lwd=lwd )
plot.by.aa( rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.3), c('unal', 'Chelicerata', 'Arachnida', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa'), g.code, lwd=lwd )

par(mfrow=c(3,1))
lwd=3
plot.by.aa( rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.3), c('unal', 'Eutheria', 'Chelicerata', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa'), g.code, lwd=lwd )
plot.by.aa( rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.3), c('unal', 'Aves', 'Chelicerata', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa'), g.code, lwd=lwd )
plot.by.aa( rbind('unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.3), c('unal', 'Mammalia', 'Chelicerata', 'Aconoidasida', 'Nematoda', 'Cnidaria', 'Lophotrochozoa'), g.code, lwd=lwd )


### We should in fact be able to fit the data as a linear model:
### The frequencies observed in the unaligned sequences may represent
### a linear combination of the ones represented in the aligned
### data set.
###
### This is somewhat sketchy though as many of these patterns
### are very similar; it's also sketchy as ideally we would
### include separate profiles for mitochondrial and non-mitochondrial
### sequences. And, so on.
### It also would seem likely that there may be many unique combinations
### that are possible.
### The selection of which codon usage profiles to include is likely to be
### important.

### Fitting the model.
### If Y is an observed set of codon frequencies in the unaligned data set
### and
### X is a set of observed codon frequencies in a set of taxons represented
### as a matrix with 1 row (n) for each codon and 1 colum (m) for each taxons
###
### then
### Y = X %*% b
###
### where b is single column matrix with m rows (1 for each taxon and column in X)
###
### from this equation we should be able to get a least squares estimate using:
###
### b = (t(X)X)^-1 %*% (t(X)Y)
###
### which can be written as:
###
### b <- solve(crossprod(X)) %*% crossprod(X,Y)
###
#### This may or may not give any reasonable inference.

## let us make a function for this. Let us simply call it 'unmix'

unmix <- function(X, Y, backSolve=FALSE){
    ## function explained amongst others at:
    ## https://rstudio-pubs-static.s3.amazonaws.com/163134_83d9fe8b2d834d18830c11a4336b7083.html#solving-system-of-equations
    if(!backSolve)
        return(solve(crossprod(X)) %*% crossprod(X,Y))
    ## However, solve() is supposed to be unstable, and the above
    ## suggests use of backsolve rather than solve doing:
    QR <- qr(X)
    Q <- qr.Q(QR)
    R <- qr.R(QR)
    backsolve( R, crossprod(Q,Y))
    ## but I haven't read up on what that does or tried it out. 
    ## Note that there is also
    ## qr.solve( )
    ## which may do all of these things.
}

## this is based on nt.codon.usage.3
o <- order(nt.codon.cor.2.3, decreasing=TRUE)
taxons <- rownames(nt.codon.usage.3)[o]

## that is a very bad list to use as it is very overlapping. It would be better to remove very similar
## profiles and to sample non-similar profiles as well. For example

plot(nt.codon.usage.3['Chelicerata',], nt.codon.usage.3['Arachnida',] )
## is almost identical;
## note that the actual numbers also indicate that most (4/5) of the 'Chelicerates' are arachnids
sum(nt.codon.usage.3['Chelicerata',]) / sum(nt.codon.usage.3['Arachnida',])

chel.b <- grepl("Chelicerata", nt.codon.u.tax )
arach.b <- grepl("Arachnida", nt.codon.u.tax )

sum(chel.b) ## 1404
sum(arach.b) ## 1180
sum(chel.b & arach.b) ## 1180 !!
## Indeed the arachnids are a subset of the chelicerates; not the other way round.
## and unfortunately most of the chelicerates here are arachnids; not the other
## way around.

## But note that all the matches to non-arachnid chelicerates appear to be to marine
## organisms. But there are only 5 such species included; two of them being different
## types of horse-shoe crabs.

## However, matches to arachnids are primarily to mites and ticks. And aquatic mites (even
## marine ones) do seem to exist.

## make an index..
mk.taxon.index <- function(){
    all.tax <- unique(nt.codon.u.tax)  ## 11518 entries
    all.tax <- strsplit(all.tax, ",")
    all.tax.nm <- unique(unlist(all.tax))
    all.tax.i <- vector(mode='list', length=length(all.tax.nm))
    names(all.tax.i) <- all.tax.nm
    for(i in 1:length(all.tax)){
        for(j in 1:length(all.tax[[i]]))
            all.tax.i[[all.tax[[i]][[j]]]] <- cbind( all.tax.i[[all.tax[[i]][[j]]]], c(i,j) )
    }
    all.tax.i
}

## a very slow and ponderous function
## uses
collect.independent.taxons <- function(taxons, all.tax.i){
    all.tax <- unique(nt.codon.u.tax)  ## 11518 entries
    all.tax <- strsplit(all.tax, ",")
    observed.tax <- c()
    collected <- c()
    for(tax in taxons){
        ## check if the taxon has been observed
        if( tax %in% observed.tax )
            next
        ## obtain all terms that are members of the taxon
        tax.members <- unique( unlist( apply( all.tax.i[[tax]], 2, function(x){
            all.tax[[x[1]]][ 1:x[2] ]
        })))
        if( !(any(tax.members %in% collected )) ){
            collected <- c(collected, tax)
            observed.tax <- c(observed.tax, tax.members )
        }
    }
    collected
}

all.tax.i <- mk.taxon.index()

sel.taxons <- collect.independent.taxons( taxons, all.tax.i )
## that gives me 33 taxons.

## we can than more or less go ahead and see if we can get a solution from these.
components <- unmix( t( nt.codon.usage.3[ sel.taxons, ] / rowSums(nt.codon.usage.3[sel.taxons,]) ), codon.usage.f2$unal[,'r_1'] / sum(codon.usage.f2$unal[,'r_1']) )

components.2 <- unmix( t( nt.codon.usage.3[ sel.taxons, ] / rowSums(nt.codon.usage.3[sel.taxons,]) ), codon.usage.f2$unal[,'r_1'] / sum(codon.usage.f2$unal[,'r_1']), backSolve=TRUE )
names(components.2) <- sel.taxons

## these look the same.. lets see what we get from that:
fitted.1 <- t( nt.codon.usage.3[ sel.taxons, ] / rowSums(nt.codon.usage.3[sel.taxons,]) ) %*% components
plot(codon.usage.f2$unal[,'r_1'], fitted.1)
## and that is a reasonable fit.
## however, we need a non-negative matrix factorisation; That is no negative components should be allowed.
## sounds like non-negative matrix factorisation, but the packages that I've found for that do different
## things (finds the a set of components in a single matrix of values)

## Non-Negative Least Squares
## Is the term that I should have been searching for. This will allow spectral unmixing. There are many implementations:
## 1. hsdar: a general package for handling hyperspectral imaging data. Not really suitable due to being
##    specialised for such data (many specific classes for meta-data and so on).
## 2. nnls: provides a very simple function that performs the operation. Very general, but maybe not so flexible.
##    does not handle weights natively, but can be simulated by multiplying the X and y by the square root of the
##    weights.
## 3. cvxr: allows weighted non-negative least squares directly (probably among other things)
## 4. nls: an R built in. Can do non-negative least squares with the 'port' algorithm with lower bounds all 0. eg:
##     zeros <- numeric(ncol(X))
##     nls(Y ~ X %*% b, start = list(b = zeros), weights = w, lower = zeros, alg = "port")
##    (though apparently not a good idea to use the zeros as a start condition?
## 5. glmnet: through setting appropriate values of alpha and gamma. Maybe a very general function.
## 6. bvls: bounded-variable least squares; where we can set the range for each value of B
##    eg: bvls(x, y, bl = rep(0, p), bu = rep(Inf, p))
##    Though we would probably set the upper bound to be 1. 
##
## options 2-4 from discussion : https://stackoverflow.com/questions/47888996/weighted-nonnegative-least-squares-in-r
## option 5 from: https://www.r-bloggers.com/2019/11/non-negative-least-squares/

## lets try with nnls, bvls, glmnet.

## use the nnls terminology
## where we solve for x
## A x = b
## b is our observed unknown data
## A is a matrix of our known profiles
require(nnls)
A <- t( nt.codon.usage.3[ sel.taxons, ] / rowSums(nt.codon.usage.3[sel.taxons,]) )
b <- codon.usage.f2$unal[,'r_1'] / sum(codon.usage.f2$unal[,'r_1'])

x.1 <- nnls( A, b )
sum(x.1$x) ## 0.9871355
plot(b, x.1$fitted)  ## very nice
par(mar=c(9.1, 4.1, 4.1, 4.1))
barplot( x.1$x, names.arg=sel.taxons, las=2 )

## is there a dependance on the order of the rows?
o <- sample(1:ncol(A))
x.1.1 <- nnls(A[,o], b)
par(mfrow=c(2,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))
##plot(b, x.1.1$fitted)  ## very nice
par(mar=c(9.1, 4.1, 4.1, 2.1))
barplot( x.1$x, names.arg=sel.taxons, las=2 )
barplot( x.1.1$x, names.arg=sel.taxons[o], las=2 )
## the result is very stable.

require(bvls)
p <- ncol(A)
x.2 <- bvls(A, b, bl=rep(0, p), bu=rep(Inf, p))
## this gives exactly the same set of values as nnls;
## suggesting that they do exactly the same thing. So not
## very useful.

require(glmnet)
x.3 <- glmnet( A, b, lambda=0, lower.limits=0, intercept=FALSE )

tmp.tax <- c('unal', 'fitted', 'Chelicerata', 'Aconoidasida', 'Amoebozoa', 'Viruses', 'Retaria')
tmp.tax <- c('unal', 'fitted', 'Chelicerata') ## , 'Aconoidasida', 'Amoebozoa', 'Viruses', 'Retaria')
tmp.cu <- list()
tmp.cu[[1]] <- rbind('fitted'=x.1$fitted[,1], 'unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.1)
tmp.cu[[2]] <- rbind('fitted'=x.1$fitted[,1], 'unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.2)
tmp.cu[[3]] <- rbind('fitted'=x.1$fitted[,1], 'unal'=codon.usage.f1$unal[,'r_1'], nt.codon.usage.3)

par(mfrow=c(3,1))
lwd=3
plot.by.aa( tmp.cu[[1]], tmp.tax, g.code, lwd=lwd )
plot.by.aa( tmp.cu[[2]], tmp.tax, g.code, lwd=lwd )
plot.by.aa( tmp.cu[[3]], tmp.tax, g.code, lwd=lwd )

cor(x.1$fitted, codon.usage.f1$unal[,'r_1']) ## 0.9349384
max( cor(codon.usage.f1$unal[,'r_1'], A) ) ## 0.91
max( cor(x.1$fitted, A) ) ## 0.968
## all are:
##      Chelicerata  Apocrita Aconoidasida  Cnidaria  Nematoda Amoebozoa   Viruses
## [1,]   0.9681523 0.9577722    0.9224155 0.9405958 0.9242387 0.9065864 0.9052604
##        Retaria unclassified entries  Mollusca Brachiopoda Clitellata Gunneridae
## [1,] 0.8784325            0.8357271 0.8207032    0.784206   0.788691  0.7759091
##      Proteobacteria Terrabacteria group Ascomycota PVC group Amphinomidae
## [1,]      0.7286934           0.6892221  0.6831694 0.6033674    0.5875348
##       Cyprinus Phyllodocida Euarchontoglires Canalipalpata uncultured eukaryote
## [1,] 0.5663977    0.5613251        0.5047191     0.4844689            0.4710952
##      Brachycera    Sauria Eupercaria   Discoba Ovalentaria Amphioxiformes
## [1,]  0.4550113 0.4457928  0.4144882 0.4160806   0.3850702      0.3272626
##      uncultured bacterium Scolecida    Clupea Stramenopiles
## [1,]            0.2864713 0.3115077 0.2444875     0.2012959

## lets do the same with al.codon.usage as well
tmp.1 <- do.call( rbind, al.codon.usage.1 )
tmp.2 <- do.call( rbind, al.codon.usage.2 )
tmp.3 <- do.call( rbind, al.codon.usage.3 )

A <- t( tmp.1[sel.taxons,] / rowSums(tmp.1[sel.taxons,]) )
b <- codon.usage.f2$unal[,'r_1'] / sum(codon.usage.f2$unal[,'r_1'])

al.x.1 <- nnls( A, b )
sum(al.x.1$x) ## 0.9740755
plot(b, al.x.1$fitted)  ##
par(mar=c(9.1, 4.1, 4.1, 4.1))
barplot( al.x.1$x, names.arg=sel.taxons, las=2 )

A <- t( tmp.2[sel.taxons,] / rowSums(tmp.2[sel.taxons,]) )
A <- A[ , !is.na(A[1,]) ]
al.x.2 <- nnls( A, b )
sum(al.x.2$x) ## 0.9854219
plot(b, al.x.1$fitted)  ##
par(mar=c(9.1, 4.1, 4.1, 4.1))
barplot( al.x.2$x, names.arg=colnames(A), las=2 )

A <- t( tmp.3[sel.taxons,] / rowSums(tmp.3[sel.taxons,]) )
A <- A[ , !is.na(A[1,]) ]
al.x.3 <- nnls( A, b )
sum(al.x.3$x) ## 0.9894078
plot(b, al.x.3$fitted)  ##
par(mar=c(9.1, 4.1, 4.1, 4.1))
barplot( al.x.3$x, names.arg=colnames(A), las=2 )

## using the query sequences rather than the database ones
## gives very variable results that are inconsistent with
## the database sequences. This suggests that we cannot
## really trust these results very much.
##
## The only other way to check this would be to bootstrap the
## nt results. Don't think that we would get much information
## either way. 
