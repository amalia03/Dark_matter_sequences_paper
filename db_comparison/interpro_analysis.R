##Protein analysis
source("functions_pep.R")
library(VennDiagram)
library(RColorBrewer)
library(gridExtra)

#Find file directory and make a variable that distinguishes between types of files. 
neg.u.f <- list.files("unmapped/orf_neg", full.names=TRUE)
neg.m.f <- list.files("mapped/orf_neg", full.names=TRUE)

type.re <- c(fa="fa$", pfam="Pfam.tsv", panther="panther.tsv", gene3d="Gene3D.tsv", superfamily="SUPERFAMILY.tsv", prints="PRINTS.tsv", prosite="ProSiteProfiles.tsv")
neg.u.f <- lapply( type.re, grep, x=neg.u.f, value=TRUE )
neg.m.f <- lapply( type.re, grep, x=neg.m.f, value=TRUE )

## Import data, exclude the first one which contains the fasta info. 
neg.u.data <- lapply( neg.u.f[-1], read.tsv.tables )
#pos.u.data <- lapply( pos.u.f[-1], read.tsv.tables )
neg.m.data <- lapply( neg.m.f[-1], read.tsv.tables )

# And also include the fasta column 
## these are the amino acid lengths of the sequences that were
## subjected to the analyes
neg.u.fa <- read.fa.l( neg.u.f[[1]] )
#pos.u.l <- read.fa.l( pos.u.f[[1]] )
neg.m.fa <- read.fa.l( neg.m.f[[1]] )

#Get the top quality matches (in terms of evalue) per query: 
neg.u.top <- top.function(neg.u.data, e.thresh=-5)
neg.m.top <- top.function(neg.m.data, e.thresh=-5)

#But given that ProSite works somewhat differently with its scoring (it is just the reverse so 9 would mean 1 in a billion (recommended threshold 8.5)), I had to do that one separately and combine it later (because anyway the Prosite entry list in the filtered list is empty).
prosite.m <- top.function.lite(neg.m.data[[6]], ecomp="greater", e.thresh=0.95, decrease=T)
prosite.u <- top.function.lite(neg.u.data[[6]], ecomp="greater", e.thresh=0.95, decrease=T)
neg.u.top[[6]] <- prosite.u
neg.m.top[[6]] <- prosite.m

##Remove the blank values(not necessary given the nature of the study)
neg.u.top.nonull <- lapply(1:length(neg.u.top), function (x){neg.u.top[[x]][neg.u.top[[x]]$description!="-",]})
neg.m.top.nonull <- lapply(1:length(neg.m.top), function (x){neg.m.top[[x]][neg.m.top[[x]]$description!="-",]})

##but keep the undescribed for reference
neg.u.blank <- lapply(1:length(neg.u.top), function (x){neg.u.top[[x]][neg.u.top[[x]]$ipro_des=="-",]})
neg.m.blank <- lapply(1:length(neg.m.top), function (x){neg.m.top[[x]][neg.m.top[[x]]$ipro_des=="-",]})

###End of data processing
###Beginning data visualization
##############
#unique GO term values for unmatched and matched list
sapply(1:length(neg.u.top), function(x){length(unique(neg.u.top[[x]]$go))})
sapply(1:length(neg.m.top), function(x){length(unique(neg.m.top[[x]]$go))})

###Plotting the top 30 counts
###
#par(mfrow=c(1,1))
#par(mar=c(4,24,4,4))

###A venn diagram function that can compare between common values ((default:protein ids) between two different databases: 
venn.d <- function(x,y,type="protid", data.l=neg.u.top){
#    jpeg(paste0(names(type.re[x+1]),"_",names(type.re[y+1]), "_shared_",type,".jpg"), width = 900, height = 900)
    grid.newpage()                    # Create new plot (works a bit differently than plot.new)
    draw.pairwise.venn(area1=length(unique(data.l[[x]][,type])),    # Draw pairwise venn diagram
                       area2=length(unique(data.l[[y]][,type])),
                       cross.area = length(unique(intersect(data.l[[x]][,type], data.l[[y]][,type]))),
                       category=c(names(type.re[x+1]),names(type.re[y+1])),
                       cat.cex=1.5, 
                       lwd=3, cex=1.3, fill=c("skyblue", "tomato"))
                                        #   dev.off()
}

#2 examples: 
#For pfam and gene3d
venn.d(1,5,type="accession", data.l=neg.u.data)
#For Pfam and superfamily
venn.d(1,6)

###Some aligned vs unaligned comparisons of the number of interpro descriptions that either shared or exclusive.  
p.db <- 6
common <- intersect(neg.u.top[[p.db]]$ipro_des, neg.m.top[[p.db]]$ipro_des)
c("matched_in_common",nrow(subset(neg.m.top[[p.db]], (neg.m.top[[p.db]]$ipro_des %in% common))))
c("matched_unique",nrow(subset(neg.m.top[[p.db]], !(neg.m.top[[p.db]]$ipro_des %in% common))))
c("unmatched_in_common",nrow(subset(neg.u.top[[p.db]], (neg.u.top[[p.db]]$ipro_des %in% common))))
c("unmatched_unique",nrow(subset(neg.u.top[[p.db]], !(neg.u.top[[p.db]]$ipro_des %in% common))))

##############################################################################################################################################################
####Get the top values for each dataset. 
top.match <- lapply(matchpoints, function(x){cbind(head(rownames(x[order(x$match, decreasing=T),]),10),
                                                   head(x[order(x$match, decreasing=T),"match"],10)
                                                   )})
top.unmatch <- lapply(matchpoints, function(x){cbind(head(rownames(x[order(x$unmatch, decreasing=T),]),10),
                                                   head(x[order(x$unmatch, decreasing=T),"unmatch"],10)
                                                   )})
###The top matching signatures were used in as a supplementary table as it provides detailed descriptions of the sequences we got. 
lapply(top.match, function(y){
cbind(sapply(y[,1], function(x){sig.index.u[sig.index.u$signature==x,2]}),y[,2])})

lapply(top.unmatch, function(y){
cbind(sapply(y[,1], function(x){sig.index.u[sig.index.u$signature==x,2]}),y[,2])})

##################################################################################################
####PLOTTING TIME

####Barplots for matches to the Interpro database
##Number of entries for each database
n.entries <- c(19000, 16000, 7000, 2000 , 2000, 1000)

###Apparently, panther has an extra subunit for each of its sections, sepaerated by a colon, if we keep the first part, then that makes for the unique group
panth.u.top <- strsplit(neg.u.top.nonull[[2]]$signature, split=":")
neg.u.top.nonull[[2]]$signature <- unlist(lapply(panth.u.top, function(x){x[1]}))

panth.m.top <- strsplit(neg.m.top.nonull[[2]]$signature, split=":")
neg.m.top.nonull[[2]]$signature <- unlist(lapply(panth.m.top, function(x){x[1]}))

###Barplots for the database comparisons 
par(mfrow=c(2,1))
par(mar=c(4,5,1,1))

###First barplot, showing the number of annotated aligned sequences
barplot(rbind(sapply(1:6, function(x){length(unique(neg.m.top.nonull[[x]]$protid))}),sapply(1:6, function(x){length(unique(neg.u.top.nonull[[x]]$protid))})),
        ylab="Number of unique counts", beside=T, cex.names=1.5, cex.lab=1.5, cex.axis=1.2, names.arg=names(type.re[-1]), yaxs="i",legend.text=(c("nt-matched", "nt-unmatch")))
box(lwd=1.5)

###Second barplot, number of unique signature IDs that aligned from the six Interpro databases. 
barplot(rbind(sapply(1:6, function(x){length(unique(neg.m.top.nonull[[x]]$signature))}),sapply(1:6, function(x){length(unique(neg.u.top.nonull[[x]]$signature))})), ylim=c(0, max(n.entries)*1.05),
        #main=" Number of unique matching Interpro signature at evalue e-10 ",
        cex.lab=1.5, cex.names=1.5, cex.axis=1.3, beside=T, col=c("gray24", "gray69"), ylab="Number of unique counts",names.arg=names(type.re[-1]))
legend("topright", c("match to nt-aligned", "match to nt-unaligned", "match to both"), col=c("gray14", "gray69", "gray89"), pch=19, cex=1.3)
rect(seq(1,16,3), 0, seq(3, 18,3), n.entries, lty=2)
box(lwd=1.5)
sapply(1:length(neg.u.top.nonull), function(x){
#    seq(1, 2*length(neg.u.top.nonull)+4, 3)
    ##Mutual values between match and unmatched
    x.offset <- seq(1, 2*length(neg.u.top.nonull)+4, 3)[x]
    rect(x.offset,0, x.offset+2,
         sum(unique(neg.m.top.nonull[[x]]$signature)%in%unique(neg.u.top.nonull[[x]]$signature)), col="gray89")
    lines(x=c(x.offset+1,x.offset+1), y=c(0, sum(unique(neg.m.top.nonull[[x]]$signature)%in%unique(neg.u.top.nonull[[x]]$signature))), lty=2)})
#############################################################################################

##PROTEIN SCATTERPLOTS
matchpoints <-lapply(2:length(type.re)-1, function(x){
   agg.col <- "signature"
#    agg.col <- "go"
    match.dl <- neg.m.top.nonull
    unmatch.dl <- neg.u.top.nonull
    all.t <- sort(unique(names(c(table(match.dl[[x]][,agg.col] ), table(unmatch.dl[[x]][,agg.col] )))))
    matchings <- data.frame(match=c(rep(0, length(all.t))),unmatch=c(rep(0, length(all.t))), row.names=all.t)
    matchings[rownames(matchings)%in%names(table(match.dl[[x]][,agg.col] )),"match"] <- table(match.dl[[x]][,agg.col] )
    matchings[rownames(matchings)%in%names(table(unmatch.dl[[x]][,agg.col] )),"unmatch"] <- table(unmatch.dl[[x]][,agg.col] )
    matchings
})

max.m <- max(sapply(1:6, function(x){max(matchpoints[[x]][,1])}))
max.u <- max(sapply(1:6, function(x){max(matchpoints[[x]][,2])}))

sig.index <- unlist(sapply(2:length(type.re)-1, function(x){c(neg.m.data[[x]][,"signature"], neg.u.data[[x]][,"signature"])}))
sig.index <- as.data.frame(cbind(sig.index,unlist(sapply(2:length(type.re)-1, function(x){c(neg.m.data[[x]][,"description"], neg.u.data[[x]][,"description"])}))), stringsAsFactors=F)
colnames(sig.index) <- c("signature", "description")
sig.index.u <- sig.index[!duplicated(sig.index$signature),]

id.l <- unlist(sapply(2:length(type.re)-1,function(i){
    unlist(sapply(1:nrow(matchpoints[[i]]),function(x){sig.index.u[sig.index.u$signature==rownames(matchpoints[[i]])[x],"description"]}))}))

############################################################################
##This was my first attempt to a protein scatterplot, below you can find the one that I actually end up using. This one I kept  because I wanted to keep the identify function as it is quite interesting.
###For the identify labels
par(mar=c(0,0,0,0), bg="white")
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE, asp=1)

hist.match <- hist(unlist(sapply(c(1:4,6), function(x){matchpoints[[x]]$match})), breaks=100, plot=F)
par(new=TRUE, plt=c(0.15, 0.8, 0.8, .95))
barplot(log(hist.match$counts+1, 10), las=2)
hist.unmatch <- hist(unlist(sapply(c(1:4,6), function(x){matchpoints[[x]]$unmatch})), breaks=100, plot=F)
par(new=TRUE, plt=c(0.8, .95, 0.15, .8))
barplot(log(hist.unmatch$counts+1, 10), horiz=T)

par(new=TRUE, plt=c(0.15, 0.8, 0.15, .8))
plot(1, type="n", ylim=c(0, max.u), xlim=c(0, max.m), xlab="nt match counts", ylab="nt unmatch counts", cex.axis=1, cex.lab=1)

sapply(c(1:4, 6), function(x){
    lines(matchpoints[[x]]$match, matchpoints[[x]]$unmatch, type="p", pch=19, cex=2, col=alpha(col.a(6)[x], .6))
})
legend("topright", legend=(names(type.re[c(2:5,7)])), pch=19, col=col.a(6)[c(1:4,6)], cex=1.4)

###Remember: when this one is executed, you can only exit it by right clicking (maybe there is another way but that is what I usually do). 
identify(
    x=unlist(sapply(2:length(type.re)-1, function(x){matchpoints[[x]]$match})),
    y=unlist(sapply(2:length(type.re)-1, function(x){matchpoints[[x]]$unmatch})),
    labels=id.l,
    plot=TRUE
)

#############################################
###Protein scatterplot, the actual one I ended up using. 
col.a <- colorRampPalette(c("hotpink2","goldenrod","peachpuff", "tomato", "firebrick3", "skyblue","skyblue4","royalblue4","mediumpurple3"))
###For the identify labels

id.pos <- cumsum(lapply(matchpoints, nrow))
matchpoints[[1]] <- cbind(matchpoints[[1]],id.l[1:nrow(matchpoints[[1]])])
matchpoints[[2]] <- cbind(matchpoints[[2]],id.l[(1+id.pos[1]):id.pos[2]])
matchpoints[[3]] <- cbind(matchpoints[[3]],id.l[(1+id.pos[2]):id.pos[3]])
matchpoints[[4]] <- cbind(matchpoints[[4]],id.l[(1+id.pos[3]):id.pos[4]])
matchpoints[[5]] <- cbind(matchpoints[[5]],id.l[(1+id.pos[4]):id.pos[5]])
matchpoints[[6]] <- cbind(matchpoints[[6]],id.l[(1+id.pos[5]):id.pos[6]])
#Remove the PRINTS database since we are not using it here..
matchpoints <- matchpoints[-5]

######################################################################
#For the ScatterPlotting
jpeg("scatterplot_all_db.jpg", width = 900, height = 900)

par(mfrow=c(1,1))
par(mar=c(0,0,0,0), bg="white")
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE, asp=1)

##Set up the border histograms with the abundance distributions. 
hist.match <- hist(unlist(sapply(c(1:length(matchpoints)), function(x){matchpoints[[x]]$match})), breaks=100, plot=F)
par(new=TRUE, plt=c(0.165, 0.825, 0.85, .95))
barplot(log(hist.match$counts+1, 10), las=2, xaxs="i")
text((length(hist.match$counts)/2), max(log(hist.match$counts+1,10))*.93, "Aligned frequency log)", cex=1.3)

hist.unmatch <- hist(unlist(sapply(c(1:length(matchpoints)), function(x){matchpoints[[x]]$unmatch})), breaks=100, plot=F)
par(new=TRUE, plt=c(0.85, .95, 0.165, .825))
barplot(log(hist.unmatch$counts+1, 10), horiz=T, xaxs="i", yaxs="i")
text(max(log(hist.unmatch$counts+1,10))*.9,length(hist.unmatch$counts)/2,"Unaligned frequency log)" , srt=270, cex=1.3)

#The main plot
par(new=TRUE, plt=c(0.15, 0.85, 0.15, .85))
plot(1, type="n", ylim=c(0, max.u), xlim=c(0, max.m), xlab="Aligned counts", ylab="Unaligned counts", cex.axis=1.3, cex.lab=1.3)

sapply(c(1:length(matchpoints)), function(x){
    lines(matchpoints[[x]]$match, matchpoints[[x]]$unmatch, type="p", pch=20+x, cex=2, col="darkgrey")
})

#Set the keywords that will be used in the coloring and legend 
##The reason why the keywords are separate is for getting the correct, distinct colors later. 
keywords <- c("Actin", "tubulin", "cysteine","trypsin", "protein kinase", "alpha/beta", "ef-hand", "transferase", "cytochrome", "p-loop")
viruses <- c("viruses", "microviridae", "spike", "capsid", "phage")
ssDNA <- c("ssdna")
ssRNA <- c("ssrna")
globin <- c("globin")

###We decided for a better ab line to be one that crosses between total matched and total unmatched. 
abline(0,length(neg.u.fa)/length(neg.m.fa), lty=2, col="darkgrey", lwd=2)

keyword.names <- c("Actin", "Tubulin", "Cysteine proteinase","Trypsin proteinase", "Protein kinase", "alpha/beta", "ef-hand", "transferase", "cytochrome","p-loop",
                   "globin", "ssDNA","ssRNA", "other viral proteins")

legend("topright", legend=c(names(type.re[c(2:5,7)]), keyword.names), pch=c(21:(21+(length(matchpoints)-1)), rep(19, length(keywords)+5)), col=c(rep("grey28", length(matchpoints)), col.a(length(keywords)), "darkgreen","cyan","deepskyblue", "lawngreen"), cex=\
1.5)

########################################################################
##########
##Finally, a Bitwise/Upsett Plot
## to make a bitwise map giving the databases in which the identifiers mapped, we can do:

neg.u.top.bf <- rep(0, length(neg.u.fa))
names(neg.u.top.bf) <- names(neg.u.fa)

neg.m.top.bf <- rep(0, length(neg.m.fa))
names(neg.m.top.bf) <- names(neg.m.fa)

db.bit <- 2^(0:5)
names(db.bit) <- type.re[-1]

for(i in 1:length(neg.u.top)){
    pid <- unique(neg.u.top[[i]]$protid)
    neg.u.top.bf[ pid ] <- bitwOr( neg.u.top.bf[ pid ], db.bit[i] )
    pid <- unique(neg.m.top[[i]]$protid)
    neg.m.top.bf[ pid ] <- bitwOr( neg.m.top.bf[ pid ], db.bit[i] )
}

neg.u.bf <- rep(0, length(neg.u.fa))
names(neg.u.bf) <- names(neg.u.fa)
neg.m.bf <- rep(0, length(neg.m.fa))
names(neg.m.bf) <- names(neg.m.fa)

for(i in 1:length(neg.u.data)){
    pid <- unique(neg.u.data[[i]]$protid)
    neg.u.bf[ pid ] <- bitwOr( neg.u.bf[ pid ], db.bit[i] )
    pid <- unique(neg.m.data[[i]]$protid)
    neg.m.bf[ pid ] <- bitwOr( neg.m.bf[ pid ], db.bit[i] )
}

u.tab <- 100 * (table( c(0:63, neg.u.bf) ) - 1) / length(neg.u.bf)
m.tab <- 100 * (table( c(0:63, neg.m.bf) ) - 1) / length(neg.m.bf)

par(mfrow=c(1,1))

bit.grid <- t(sapply(0:63, function(x){ 0 + as.logical(bitwAnd( x, 2^(0:5))) }))
##Testing
sapply(db.bit, function(x){ sum(as.logical(bitwAnd( neg.m.bf, x ))) / length(neg.m.bf) })
a <- cbind(round(m.tab,2),round(u.tab,2), bit.grid)
colSums(a[rowSums(a[,3:8])==5,])

#####
###Plotting the upsett
jpeg("protein_db_distr_ord.jpg", height = 900, width = 1260)
#jpeg("protein_db_distr_ratio.jpg", width = 900, height = 1260)

par(mar=c(0,0,0,0), bg="white")
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE, asp=1)

###Sort the values
m.u.ord <- sort(round(u.tab)+round(m.tab), decreasing=T)
m.u.ord <- as.numeric(names(m.u.ord[m.u.ord>0]))

##For the grid
par(new=TRUE, plt=c(0.1, 0.9, 0.1, .35))

plot(0, 0, type="n", ylim=c(0,length(db.bit)), xlim=c(0, length(m.u.ord)+1), xlab="", ylab="", axes=F,xaxs="i", yaxs="i")


u.perc <- sapply(1:length(db.bit), function(x){
    round((length(unique(neg.u.data[[x]]$protid))/length(neg.u.fa))*100)})

m.perc <- sapply(1:length(db.bit), function(x){
    round((length(unique(neg.m.data[[x]]$protid))/length(neg.m.fa))*100)})

sapply(0:(length(m.u.ord)-1),function(j){
    sapply(0:(length(db.bit)-1),function(i){
#        rect(i,j,i+1,j+1, col=ifelse(bit.grid[m.u.ord[j+1]+1,i+1]==1,"grey","white"))
        rect(j,i,j+1,i+1, col=ifelse(bit.grid[m.u.ord[j+1]+1,i+1]==1,"grey","white"))
    })
})

text(.5,seq(0.5,length(db.bit)),paste0("m= ",m.perc,"%","\n","u= ", u.perc,"%"))

text(length(m.u.ord)+.5,seq(0.5,length(db.bit)),
     labels=paste0(names(type.re[-1])),
     cex=1.3)

##For the barplot
par(new=TRUE, plt=c(0.1, 0.9, 0.36, .9))
plot(0, 0, type="n", ylim=c(0,100-min(c(m.tab[1],u.tab[1]))), xlim=c(0, length(m.u.ord)+1), xlab="", axes=F, main="", yaxs="i", xaxs="i", ylab="Proportion (%)")

axis(2,at=seq(0,80,10),labels=seq(0,80,10), col="darkgrey", cex=.85, las=2) 
abline(v=1:length(m.u.ord),h=seq(10,80,10), lty=2, col="lightgrey")
abline(v=1, lty=2, col="darkgrey")

rect(0.1,0,.5,100-u.tab[m.u.ord[1]+1], col="lightblue4")
rect(.5,0,.9,100-m.tab[m.u.ord[1]+1], col="lightblue")

sapply(1:(length(m.u.ord)-1),function(j){
    rect(j+0.1,0,j+.5,u.tab[m.u.ord[j+1]+1], col="lightblue4")
    })

box(col="darkgrey")
sapply(1:(length(m.u.ord)-1),function(j){
    rect(j+.5,0,j+.9,m.tab[m.u.ord[j+1]+1], col="lightblue")
    })
legend("topright", legend=c("nt unmatch", "nt match"), col=c("lightblue4", "lightblue"), lty=1, lwd=3, bg="white", cex=1.5)

dev.off()

################
###################
##Supplement: Getting the actin matches. 

#sapply(1:length(db.bit), function(x){nrow((neg.m.top[[x]][neg.m.top[[x]]$ipro_des=="Actin family",]))})
actin <- unique(c(unlist(sapply(1:length(db.bit), function(x){neg.m.top[[x]][neg.m.top[[x]]$ipro_des=="Actin family","protid"]}))),
                unlist(sapply(1:length(db.bit), function(x){neg.u.top[[x]][neg.u.top[[x]]$ipro_des=="Actin family","protid"]})))

write.table(actin, "actin_matches.tsv", sep="\t", row.names=F, col.names=F, quote=F)
