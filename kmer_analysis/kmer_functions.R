#!/usr/bin/ Rscript

dyn.load("/directory/to/C_kmer_functions/R_c_plugins/nuc_kmer_count.so")

kmerify <- function(seqs,k.size){
    counts <- .Call("count_kmers", seqs, k.size)
    kmers <- .Call("ints_to_kmers", as.integer( (1:length(counts))-1 ), k.size)
    kmer.df <- data.frame(kmers, counts)
    kmer.df
}

bin.range <- function(b, range=a.u.range){
    a <- sapply(1:length(range), function(x){
        bin <- rep(0,length(b))
        bin[which(b>=range[x]& b<=range[x+1])] <- x
        bin
    })
    rowSums(a)
}

#Input here must be a kmer df in the format outputed by the kmerify function. 

dens.pixel.plot <- function(kmers.a, kmers.b, pixlevel=70, x.name="", y.name="", mt="", col1="royalblue4", col2="tomato"){
##In case there is a graph loop, it is better to reset the parameters
    old.par <- par(mar = c(0, 0, 0, 0))
    par(old.par)

    a.u.range <- seq(min(c(kmers.b,kmers.a)),max(c(kmers.b, kmers.a)),(max(c(kmers.b, kmers.a))-min(c(kmers.b,kmers.a)))/pixlevel)

    ab <- round(log(table(bin.range(kmers.a, range=a.u.range),bin.range(kmers.b, range=a.u.range)),10)+1,1)

    u.ab <- sort(unique(ab))
    u.ab[which(!is.finite(u.ab))] <- 0

    pal2<- colorRampPalette(c(col1,col2))(round(max(u.ab), digits=0))

   plot(min(a.u.range): max(a.u.range),min(a.u.range): max(a.u.range), type="n" , ylim=c(min(a.u.range), max(a.u.range)), xlab=x.name, ylab=y.name, xlim=c(min(a.u.range), max(a.u.range)), main=mt)
    sapply(1:length(a.u.range),function(x){
        sapply(1:length(a.u.range),function(y){
            rect(a.u.range[y],
                 a.u.range[x],
                 a.u.range[y+1],
                 a.u.range[x+1],border="NA",
                 col=pal2[ab[as.numeric(rownames(ab))==x,as.numeric(colnames(ab))==y]])
        })
    })
    ##For the legend

    par(new=TRUE, plt=c(0.15, 0.3, 0.75, 0.88))
    plot(0:1, 0:1,type="n", axes=F, ylab="", xlab="", xaxs="i",ylim=c(0,length(pal2)*1.6))
    rect(0.1,
         1:length(pal2),
         0.5,
         2:(length(pal2)+1),
         col=pal2[seq(1, length(pal2),1)],
         border=pal2
         )        
           text(x=0.8, y=c(1.5, 1+length(pal2)/2, length(pal2)+0.5),
         labels=c(u.ab[2], u.ab[round(length(u.ab)/2, digits=0)], u.ab[length(u.ab)]))

    text(x=0.5, y=length(pal2)*1.5, labels="Density (log 10)", cex=0.80)
}

symmetrify.dens <- function(seqs, k.size, type, maint="", col1="deepskyblue", col2="deeppink3", pix=100){
    counts <- .Call("count_kmers", seqs, k.size)
    kmers <- .Call("ints_to_kmers", as.integer( (1:length(counts))-1 ), k.size)
    kmer.df <- data.frame(kmers, counts)
    r.counts <- .Call("rc_kmer_ints", as.integer( (1:length(counts))-1 ), k.size)
    rev.comp <- kmer.df[r.counts+1,]
    kmer.l <- round(log(kmer.df$counts,10), digits=2)
    rc.l <- round(log(rev.comp$counts,10), digits=2)
    dens.pixel.plot(kmer.l, rc.l, pixlevel=pix, y.name="kmer count(log10)", x.name="reverse complement kmer (log10)", mt=maint, col1=col1, col2=col2)
}


