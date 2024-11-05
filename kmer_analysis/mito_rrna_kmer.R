###kmer analysis using lmjs kmer script (more efficient)
library("scales")
dyn.load("/home/ama/R/R_c_plugins/nuc_kmer_count.so")

al.fa <- readLines("all_aligned_300.fa")
unal.fa <- readLines("all_unaligned_300.fa")

mitos <- readLines("mito_rrna.fa")
no.mitos <- readLines("no_mito_rrna.fa")

####Last one is an empty string so I am removing that one
al.fa <- al.fa[-length(al.fa)]
unal.fa <- unal.fa[-length(unal.fa)]

mitos <- mitos[-length(mitos)]
no.mitos <- no.mitos[-length(no.mitos)]

al.ids <- al.fa[seq(1,length(al.fa)-1, 2)]
unal.ids <- unal.fa[seq(1,length(unal.fa)-1, 2)]
al.seqs <- al.fa[seq(2,length(al.fa), 2)]
unal.seqs <- unal.fa[seq(2,length(unal.fa), 2)]

mito.seqs <- mitos[seq(2,length(mitos),2)]
no.mito.seqs <- no.mitos[seq(2,length(no.mitos),2)]

k.s <- 7L
pal2<- colorRampPalette(c("navy","magenta"))

###Plotting aligned vs unaligned
#sapply(5:10,function(k.s){
#jpeg(paste0("mito_rrna_vs_unal_kmer_", k.s, ".jpg"), width = 960, height = 480)
###jpeg(paste0("kmer_miseq_al_vs_unal_kmer_", k.s, "_logged.jpg"), width = 960, height = 480)
    par(mfrow=c(1,1))
    write(paste("Starting analysis with kmer size =", k.s, "at", Sys.time()), stdout())
    kmerify <- function(seqs,k.size,r.c=0){
        counts <- .Call("count_kmers", seqs, k.size)
        kmers <- .Call("ints_to_kmers", as.integer( (1:length(counts))-1 ), k.size)
        kmer.df <- data.frame(kmers, counts)
        kmer.df
    }

al.kmr <- kmerify(al.seqs,k.s, r.c=0)
    unal.kmr <- kmerify(unal.seqs,k.s, r.c=0)
    mito.kmr <- kmerify(mito.seqs,k.s,r.c=0)
    no.mito.kmr <- kmerify(no.mito.seqs,k.s,r.c=0)
    mito.pct <- round((mito.kmr$counts/al.kmr$counts)*100,0)
    col.pal<- data.frame("pct"= seq(0, 100 ,1),"col.index"=as.character(pal2(101)))
    al.kmr$col <- sapply(1:length(mito.pct), function(x){
        col.pal[col.pal$pct==mito.pct[x], "col.index"]
    })
    plot(log(al.kmr$counts,10), log(unal.kmr$counts,10),  ylab="Unaligned reads (log10)", xlab="Aligned reads (log 10)", main=paste("Aligned Vs Unaligned miseq reads at kmer =", k.s), pch=19,col=alpha(al.kmr$col,0.5))
    legend("topleft",c("kmer match to mito and/or rRNA(%)","kmer not match to mito or rRNA(%)"),text.col=c("magenta","royalblue"))
#dev.off()
#})

lin.regr <- function(data.a, data.b){
    r2 <- summary(lm(log(data.a,10)~log(data.b,10)))
    r2 <- round(r2$r.squared,digits=3)
    r2
}

###Plotting before and after...
#sapply(5:10,function(k.s){
#jpeg(paste0("mito_rrna_vs_unal_kmer_", k.s, ".jpg"), width = 960, height = 480)
###jpeg(paste0("kmer_miseq_al_vs_unal_kmer_", k.s, "_logged.jpg"), width = 960, height = 480)
    par(mfrow=c(1,2))
    write(paste("Starting analysis with kmer size =", k.s, "at", Sys.time()), stdout())
    al.kmr <- kmerify(al.seqs,k.s, r.c=0)
    unal.kmr <- kmerify(unal.seqs,k.s, r.c=0)
    mito.kmr <- kmerify(mito.seqs,k.s,r.c=0)
    no.mito.kmr <- kmerify(no.mito.seqs,k.s,r.c=0)
    mito.pct <- round((mito.kmr$counts/al.kmr$counts)*100,0)
    col.pal<- data.frame("pct"= seq(0, 100 ,1),"col.index"=as.character(pal2(101)))
    al.kmr$col <- sapply(1:length(mito.pct), function(x){
        col.pal[col.pal$pct==mito.pct[x], "col.index"]
    })

plot(log(al.kmr$counts,10), log(unal.kmr$counts,10),  ylab="Unaligned reads (log10)", xlab="Aligned reads (log 10)", main=paste("Aligned Vs Unaligned miseq reads at kmer =", k.s), pch=19,col=alpha(al.kmr$col,0.1), cex.axis=1.4, cex.lab=1.4)

#abline(0,1, col="darkgrey", lty=2, lwd=2)
abline(0,sum(log(unal.kmr$counts,10))/ sum(log(al.kmr$counts,10)),col="darkgrey", lty=2, lwd=2)

#with(par(),text(usr[2],usr[4],labels=paste("R2=", lin.regr(al.kmr$counts, unal.kmr$counts)), cex=1, font=2, adj=c(1.2,1.4)))
p.v <- round(as.numeric(cor.test(al.kmr$counts, unal.kmr$counts)[3]),2)
paste("p=", ifelse(p.v==0, ">0.001", p.v))), cex=1.8, bty="n", border="white")
legend("topright",c(paste("r=",round(cor(al.kmr$counts, unal.kmr$counts),2)),
                    paste("p=", ifelse(p.v==0, ">0.001", p.v))),
                    cex=1.7, bty="n", border="white")

plot(log(no.mito.kmr$counts,10), log(unal.kmr$counts,10),  ylab="Unaligned reads (log10)", xlab="Aligned reads, no mito (log 10)", main=paste("Aligned Vs Unaligned miseq reads at kmer =", k.s), pch=19,col=alpha("navy",0.1), cex.axis=1.4, cex.lab=1.4)
#abline(0,1, col="darkgrey", lty=2, lwd=2)

abline(0,sum(log(unal.kmr$counts,10))/ sum(log(no.mito.kmr$counts,10)),col="darkgrey", lty=2, lwd=2)
with(par(),text(usr[2],usr[4],labels=paste("R2=", lin.regr(no.mito.kmr$counts, unal.kmr$counts)), cex=1, font=2, adj=c(1.2,1.4)))

p.v <- round(as.numeric(cor.test(no.mito.kmr$counts, unal.kmr$counts)[3]),2)
paste("p=", ifelse(p.v==0, ">0.001", p.v))), cex=1.8, bty="n", border="white")
legend("topright",c(paste("r=",round(cor(no.mito.kmr$counts, unal.kmr$counts),2)),
                    paste("p=", ifelse(p.v==0, ">0.001", p.v))),
                    cex=1.7, bty="n", border="white")

###Trying some other versions of the graph
al.kmr.nt <- al.kmr[-which(al.kmr$counts==max(al.kmr$counts)),]
unal.kmr.nt <- unal.kmr[-which(al.kmr$counts==max(al.kmr$counts)),]
nomito.kmr.nt <- no.mito.kmr[-which(al.kmr$counts==max(al.kmr$counts)),]

plot(log(al.kmr.nt$counts,10), log(unal.kmr.nt$counts,10),  ylab="Unaligned reads (log10)", xlab="Aligned reads (log 10)", main=paste("Aligned Vs Unaligned miseq reads at kmer =", k.s), pch=19,col=alpha(al.kmr$col,0.1))

abline(0,1, col="darkgrey", lty=2, lwd=2)
with(par(),text(usr[2],usr[4],labels=paste("R2=", lin.regr(al.kmr.nt$counts, unal.kmr.nt$counts)), cex=1, font=2, adj=c(1.2,1.4)))

plot(log(nomito.kmr.nt$counts,10), log(unal.kmr.nt$counts,10),  ylab="Unaligned reads (log10)", xlab="Aligned reads (log 10)", main=paste("Aligned Vs Unaligned miseq reads at kmer =", k.s), pch=19,col=alpha("navy",0.1))

plot(nomito.kmr.nt$counts, unal.kmr.nt$counts,  ylab="Unaligned reads (log10)", xlab="Aligned reads (log 10)", main=paste("Aligned Vs Unaligned miseq reads at kmer =", k.s), pch=19,col=alpha("navy",0.1))

abline(lm(log(unal.kmr.nt$counts,10)~log(nomito.kmr.nt$counts,10)), col="darkgrey", lty=2, lwd=2)
abline(0,1, col="darkgrey", lty=2, lwd=2)

with(par(),text(usr[2],usr[4],labels=paste("R2=", lin.regr(nomito.kmr.nt$counts, unal.kmr.nt$counts)), cex=1, font=2, adj=c(1.2,1.4)))
