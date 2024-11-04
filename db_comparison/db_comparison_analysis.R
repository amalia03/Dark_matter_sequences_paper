#load("../.Rdata")
library(gplots)
#Import fasta file used for mutual comparison between groups

#First, a 300k subset of the original fasta file.  
fasta.file <- readLines("assembled_300K.fa")
fa.ids <- fasta.file[seq(1,length(fasta.file), 2)]
fa.ids <- unlist(strsplit(fa.ids," "))[seq(1,length(fa.ids)*2,2)]
fa.ids <- unlist(strsplit(fa.ids,">"))[seq(2,length(fa.ids)*2,2)]


#Import nr data
nr.bl <- read.delim("blastx_vs_300bp.bl", header=F)

nr.match <- unique(nr.bl[,1])
nr.match <- as.character(nr.match)

nr.match <- nr.match[nr.match%in%fa.ids]
length(nr.match)

###Sequences aligned in Rfam using Infernal
inf.ids <- readLines("/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/protein_discovery/infernal/inf_q.tsv")
inf.ids <- unique(inf.ids)
inf.match <- inf.ids[inf.ids%in%fa.ids]

###Sequences that were NOT identified by MG-RAST (not sure why but here I am importing a negative of a negative, i e what was not found is not included)
mg.unmatch <- readLines("/DATA/MarineGenomics/ama/data/miseq_reassembly/orf_analysis/miseq_samples/protein_discovery/mg_rast/mgrast_dark_matter.fa")

mg.ids <- mg.unmatch[seq(1,length(mg.unmatch), 2)]
mg.ids <- gsub(">(\\S*)","\\1", mg.ids)
mg.ids <- gsub("_.*","", mg.ids)

mg.match <- fa.ids[!(fa.ids%in%mg.ids)]

###Uniprot matches
uni.match <- read.delim("../../uniprot/all_uniprot_annotated2.bl", sep="\t", quote="", header=F, stringsAsFactors=F)

unipr.match <- unique(uni.match$V1)

length(unipr.match)
upro.match <- unipr.match[unipr.match%in%fa.ids]

#Sequences that aligned to some database in Interpro
ipro.u.match <- unique(unlist(lapply(1:length(db.bit), function(x){neg.u.data[[x]]$protid})))
ipro.m.match <- unique(unlist(lapply(1:length(db.bit), function(x){neg.m.data[[x]]$protid})))
ipro.match <- unique(c(ipro.m.match,ipro.u.match))
ipro.match <- ipro.match[ipro.match%in%fa.ids]
length(ipro.match)

##nt ids, the ones that matched
nt.match <- unique(m.len$query)
nt.match <- nt.match[nt.match%in%fa.ids]


######

#####
fastas.match <- rep(0, length(fa.ids))
names(fastas.match) <- fa.ids

db.combo <- list(nr.match, nt.match, ipro.match, upro.match, inf.match, mg.match)
dbs.bit <- 2^(0:(length(db.combo)-1))

names(dbs.bit) <- c("nr", "nt", "interpro", "uniprot", "rfam", "mg-rast")

###One thing i dont understand here is how the ids know recognition happen here...
for (x in 1:length(dbs.bit)){
    fastas.match[db.combo[[x]]] <- bitwOr( fastas.match[db.combo[[x]]], dbs.bit[x])
    }

length(db.combo)

fast.tab <- 100 * (table( c(0:(2^length(db.combo)-1), fastas.match) ) - 1) / length(fastas.match)

#bit.grid <- t(sapply(0:7, function(x){ 0 + as.logical(bitwAnd( x, 2^(0:2))) }))

bit.grid <- t(sapply(0:(2^length(db.combo)-1), function(x){ 0 + as.logical(bitwAnd( x, 2^(0:(length(db.combo)-1))))}))

jpeg("all_db_distr.jpg", height = 900, width = 1260)
dev.off()

par(mar=c(0,0,0,0), bg="white")
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE, asp=1)
fast.ord <- sort(round(fast.tab,2), decreasing=T)
fast.ord <- as.numeric(names(fast.ord[fast.ord>0]))
par(new=TRUE, plt=c(0.1, 0.9, 0.1, .35))

plot(0, 0, type="n", ylim=c(0,length(dbs.bit)), xlim=c(0, length(fast.ord)+2), xlab="", ylab="", axes=F,xaxs="i", yaxs="i")                                                                                                             \

bit.grid.ord <- bit.grid[fast.ord+1,]

sapply(0:(nrow(bit.grid.ord)-1), function(j){
    if(sum(bit.grid.ord[j+1,])==1){
        sapply(0:(ncol(bit.grid.ord)-1), function(i){
            rect(j,i,j+1,i+1, col=ifelse(bit.grid.ord[j+1,i+1]==1,"lightblue4","white"))
            })
#            rect(i,j,i+1,j+1, col=ifelse(a[i,j]==1,"red","white"))
    }else{
                sapply(0:(ncol(bit.grid.ord)-1), function(i){
        rect(j,i,j+1,i+1, col=ifelse(bit.grid[fast.ord[j+1]+1,i+1]==1,"grey","white"))})
                }})


#rect(0,0, 1, length(dbs.bit), col=alpha("khaki",.4))
rect(0,0, 1, length(dbs.bit), col=("white"))
rect(0,0, 1, (length(unique(unlist(db.combo)))/length(fa.ids))*length(dbs.bit), srt=90, col=("lightblue4"))

axis(2, at=c(seq(0,length(dbs.bit),length(dbs.bit)/5)), labels=seq(0,100,20), outer="Hello", las=2)

db.perc <- sapply(1:length(db.combo), function(x){round(length(db.combo[[x]])/length(fa.ids)*100)})

par(new=TRUE, plt=c(0.9, 1, 0.1, .35))
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE, asp=1, ylim=c(0,length(dbs.bit)))

#text(.1,seq(par()$usr[3], par()$usr[4], par()$usr[4]/length(dbs.bit)),
text(.2, seq(0.5,length(dbs.bit)+.5,1),
          paste0(names(dbs.bit),"= ", db.perc,"%"), cex=.85)

##For the barplot
par(new=TRUE, plt=c(0.1, 0.9, 0.365, .9))
#plot(0, 0, type="n", ylim=c(0,(100-max(fast.tab))), xlim=c(0, length(fast.ord)+1), xlab="", axes=F, main="", yaxs="i", xaxs="i", ylab="Proportion (%)")
plot(0, 0, type="n", ylim=c(0,max(fast.tab[-1]*1.07)), xlim=c(0, length(fast.ord)+2), xlab="", axes=F, main="", yaxs="i", xaxs="i", ylab="Proportion (%)")

axis(2,at=seq(0,80,2),labels=seq(0,80,2), col="darkgrey", cex=.85, las=2)
abline(v=1:length(fast.ord),h=seq(0,80,1), lty=2, col="lightgrey")

#rect(0.1, 0, .9, 100- (fast.tab[1]), col="lightblue")
sapply(1:(length(fast.ord)-1),function(j){
    rect(j+0.1,0,j+.9,fast.tab[fast.ord[j+1]+1], col="lightblue4")
})

box(col="darkgrey")

dev.off()

###############
##To find the proportion of sequences that did not match to MG-RAST
##Remember the dbs.bits first, which one fits with which
dbs.bit
##For MG Rast
sapply(dbs.bit, function(x){sum(fast.tab[bitwAnd(x, 0:63)>0])})

##For Interpro, Uniprot and nr combined
funct.gr <- fast.tab[bitwAnd(1, 0:63)>0 | bitwAnd(4, 0:63)>0 | bitwAnd(8, 0:63)>0]
sum(funct.gr)

##Nor functional ones
no.funct.gr <- fast.tab[bitwAnd(2, 0:63)>0 | bitwAnd(16, 0:63)>0]
sum(no.funct.gr)

no.mg.gr <- fast.tab[bitwAnd(1, 0:63)>0 | bitwAnd(2, 0:63)>0 | bitwAnd(4, 0:63)>0 |bitwAnd(8, 0:63)>0 | bitwAnd(16, 0:63)>0]
sum(no.mg.gr)

sum(funct.gr)

##For Interpro, Uniprot and nr &MGRAST combinedcombined
sum(fast.tab[bitwAnd(1, 0:63)>0 | bitwAnd(2, 0:63)>0 | bitwAnd(8, 0:63)>0]| bitwAnd(8, 0:63)>0])

##And to exclude the MG RAST
##From the protein sequences
sum(no.mg.gr[names(no.mg.gr)%in%names(fast.tab[bitwAnd(32, 0:63)==0])])
sum(funct.gr[names(funct.gr)%in%names(fast.tab[bitwAnd(32, 0:63)==0])])
sum(no.funct.gr[names(no.funct.gr)%in%names(fast.tab[bitwAnd(32, 0:63)==0])])

#################
##The uniqueness index

#a <- c(1:6)
#b <- c(2,3,5,6)
#cc <- c(3,6,9,10,11)
#d <- c(11:15)
#e <- c(1:6,50:90)
##all.l <- list(a,b,cc,d,e)
#all.t <- table(c(a,b,cc,d,e))

db.div.indx <- function(db){
    all.l <- db
    all.t <- table(c(unlist(db)))
    lapply(1:length(all.l), function(x){
        sum.a <- sum(length(all.l)/all.t[names(all.t) %in% all.l[[x]]])
        nonoff <- sum.a/(length(all.l)*length(all.l[[x]]))
        min.l <- length(all.l[[x]])/(length(all.l[[x]])*length(all.l))
        round((nonoff-min.l)*(1/(1-min.l)),2)
    })
}

###So this index is useful in providing an index of unique values relative to another group. It does not however provide information about the match size of each compared group.


db.combo <- list(nt.match, nr.match, upro.match,ipro.match, inf.match, mg.match)
names(db.combo) <- c("nt", "nr", "uniprot", "interpro","rfam", "mg")

a <- sapply(1:length(db.combo), function(x){
    sapply(1:length(db.combo), function(y){
        db <- list(db.combo[[x]],db.combo[[y]])
        db.div.indx(db)[1]
    })
})

                                        #a <- t(a)
db <- list(db.combo[[1]],db.combo[[3]])
db.div.indx(db)

mat.a <- matrix(as.numeric(a), ncol=length(db.combo), nrow=length(db.combo))
colnames(mat.a) <- names(db.combo)

rownames(mat.a) <- c(paste0("vs_", names(db.combo[1])),
                     paste0("vs_", names(db.combo[2])),
                     paste0("vs_", names(db.combo[3])),
                     paste0("vs_", names(db.combo[4])),
                     paste0("vs_", names(db.combo[5])),
                     paste0("vs_", names(db.combo[6])))


########
##########


breaks <- seq(0,max(mat.a), max(mat.a)/20)
cols <- colorRampPalette(c("sandybrown", "maroon", "lightskyblue4"))(length(breaks)-1)
#cols <- colorRampPalette(c( "sandybrown","tomato2","red"))(length(breaks)-1)
                                        #For all
par(mar=c(5,6,5,5))
plot(1:10, type="n", xlim=c(0,ncol(mat.a)), ylim=c(0, nrow(mat.a)), axes=F, ylab="",yaxs="i")
axis(2, 0.5:(nrow(mat.a)-0.5),labels=rownames(mat.a), las=2, cex=.6, tick=F)
axis(1,0.5:(ncol(mat.a)-.5), labels=colnames(mat.a),  cex=.6, tick=F)

sapply(1:nrow(mat.a), function(x){
        col.vec <- cut(mat.a[,x], breaks, labels = cols)
            rect(x-1, 0:nrow(mat.a), x, 0:nrow(mat.a)+1, col=as.character(col.vec), border=as.character(col.vec))})

sapply(1:nrow(mat.a), function(x){
    #text(x+.5,x:nrow(tax.group.h)+.5, labels=tax.group.h[x,tax.group[1:nrow(tax.group),], col="white"])
    text(x-.5, 0.5:nrow(mat.a), mat.a[,x], col="white")
    })

