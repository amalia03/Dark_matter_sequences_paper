###kmer analysis using lmjs kmer script (more efficient)
library("scales")
library("vegan")
library("MASS")
source("symmetry_functions.R")

dyn.load("/home/ama/R/R_c_plugins/nuc_kmer_count.so")

##
geans.files <- list.files("GEANS/", pattern="nless.fa$", recursive=T, full.names=TRUE)
geans.names <- sub("GEANS//(\\w+_[a-z]+)_.+", "\\1", geans.files)

l.geans <- lapply(geans.files,function(i){
    #read.table(i, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
    readLines(i)
})
names(l.geans) <- geans.names

l.geans.seqs <- lapply(1:length(l.geans), function(x){
    l.geans[[x]][seq(2,length(l.geans[[x]]), 2)]
}
)

k.s=7L
geans.kmr <- lapply(1:length(l.geans.seqs), function(x){
    kmerify(l.geans.seqs[[x]],k.s)
})
names(geans.kmr) <- geans.names


###To make an nmds
geans.kmr.df <- do.call(cbind, lapply(geans.kmr, function(x){x[,2]}))
swiss.dist <- vegdist(log(t(geans.kmr.df)), method="bray")
swiss.mds <- metaMDS(swiss.dist)
col.r <- colorRampPalette(c("tomato", "dodgerblue"))(4)
plot(swiss.mds$points, pch=16, cex=1.3, col=sapply(1:4, function(x){rep(col.r[x], 6)}))

gr.matrix <- t(sel.matrix)
if(not.all==TRUE){
    rownames(gr.matrix) <- sample.names}
else if(not.all==FALSE){
        rownames(gr.matrix) <- sample.names[1:9]}
swiss.dist <- vegdist(gr.matrix, method="bray")
swiss.mds <- isoMDS(swiss.dist)
swiss.mds
col.r <- colorRampPalette(c("tomato", "dodgerblue"))(4)

### To make a PCA plot 
#sapply(7:9, function(k.s){
    write(paste("Starting analysis with kmer size =", k.s, "at", Sys.time()), stdout())

geans.kmr <- lapply(1:length(l.geans.seqs), function(x){
    kmerify(l.geans.seqs[[x]],k.s)
})


names(geans.kmr) <- geans.names
kmr.t <- t(geans.kmr[[1]])

colnames(kmr.t) <- kmr.t[1,]
kmr.t <- as.numeric(kmr.t[2,])

2:length(geans.kmr)
for(i in 2:length(geans.kmr)){
    kmr.t <- rbind(as.numeric(t(geans.kmr[[i]])[2,]), kmr.t)
}

dim(kmr.t)
rownames(kmr.t) <- rev(geans.names)

kmr.pr <- prcomp(kmr.t, scale=T )
eigs <- kmr.pr$sdev^2

eig.values <- rbind(
    SD = sqrt(eigs),
    Proportion = eigs/sum(eigs)*100,
    Cumulative = cumsum(eigs)/sum(eigs))

#jpeg(paste0("pcas_GEANS.jpg"), width = 980, height = 900)
#par(mfrow=c(1,3), mar=c(6,6,1,1))
par(mar=c(5,6,1,1))
k.s <- 7L

col.r <- colorRampPalette(c("tomato", "dodgerblue"))(4)
plot(kmr.pr$x, col=col.r[rep(c(1,2,3,4),each=6)], pch=16, cex=3,
     xlim=c(min(kmr.pr$x[,1])*1.1,max(kmr.pr$x[,1])*1.1), xlab=paste0("PC1 ", round(eig.values[2,1],2),"%",2), ylab=paste0("PC2 ", round(eig.values[2,2], 2), "%"),
     ylim=c(min(kmr.pr$x[,2]*1.1),max(kmr.pr$x[,2]*1.1)), cex.axis=2.5, cex.lab=2.5
     )
abline(h=seq(-300,300,50),v=seq(-300,300,50), lty=2, col="darkgrey", lwd=2)
legend("bottomleft", legend=c("120", "330", "840", "ZVL"), pch=16, cex=3, col=col.r, bg="white", title="Location")
text(kmr.pr$x, labels=rep(c("T", "B"),12), pos=sample(1:4, size=nrow(kmr.pr$x), replace=T), cex=2.5)

#})

dev.off()

