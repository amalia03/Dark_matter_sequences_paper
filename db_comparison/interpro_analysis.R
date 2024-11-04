source("functions_pep.R")
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(eulerr)


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
##Total unmapped fasta sequences
#1621657

neg.u.f <- list.files("unmapped/orf_neg", full.names=TRUE)
#pos.u.f <- list.files("unmapped/orf_pos", full.names=TRUE)
neg.m.f <- list.files("mapped/orf_neg", full.names=TRUE)
type.re <- c(fa="fa$", pfam="Pfam.tsv", gene3d="Gene3D.tsv", prints="PRINTS.tsv", prosite="ProSiteProfiles.tsv", supf="SUPERFA
MILY.tsv", panther="panther.tsv")


type.re <- c(fa="fa$", pfam="Pfam.tsv", panther="panther.tsv", gene3d="Gene3D.tsv", superfamily="SUPERFAMILY.tsv", prints="PRI
NTS.tsv", prosite="ProSiteProfiles.tsv")

neg.u.f <- lapply( type.re, grep, x=neg.u.f, value=TRUE )
#pos.u.f <- lapply( type.re, grep, x=pos.u.f, value=TRUE )
neg.m.f <- lapply( type.re, grep, x=neg.m.f, value=TRUE )

## exclude the first one which contains fasta data.n
neg.u.data <- lapply( neg.u.f[-1], read.tsv.tables )
#pos.u.data <- lapply( pos.u.f[-1], read.tsv.tables )
neg.m.data <- lapply( neg.m.f[-1], read.tsv.tables )

#Get the top match of a query
neg.u.top <- top.function(neg.u.data, e.thresh=-5)
#pos.u.top <- top.function(pos.u.data, e.thresh=-10)
neg.m.top <- top.function(neg.m.data, e.thresh=-5)

#But given that ProSite works somewhat differently with its scoring (it is just the reverse so 9 would mean 1 in a billion (re
commended threshold 8.5)), I had to do that one separately and combine it later (because anyway the Prosite entry list in the
filtered list is empty).
prosite.m <- top.function.lite(neg.m.data[[6]], ecomp="greater", e.thresh=0.95, decrease=T)
prosite.u <- top.function.lite(neg.u.data[[6]], ecomp="greater", e.thresh=0.95, decrease=T)

neg.u.top[[6]] <- prosite.u
neg.m.top[[6]] <- prosite.m

##Remove the blank valeues(not necessary given the nature of the study)
#Later edit, skip this step
#neg.u.top <- lapply(1:length(neg.u.top), function (x){neg.u.top[[x]][neg.u.top[[x]]$ipro_des!="-",]})
#pos.u.top <- lapply(1:length(pos.u.top), function(x){neg.u.top[[x]][neg.u.top[[x]]$ipro_des!="-",]})
neg.u.top.nonull <- lapply(1:length(neg.u.top), function (x){neg.u.top[[x]][neg.u.top[[x]]$description!="-",]})
#pos.u.top <- lapply(1:length(pos.u.top), function(x){neg.u.top[[x]][neg.u.top[[x]]$ipro_des!="-",]})
neg.m.top.nonull <- lapply(1:length(neg.m.top), function (x){neg.m.top[[x]][neg.m.top[[x]]$description!="-",]})

##but keep the undescribed for reference
neg.u.blank <- lapply(1:length(neg.u.top), function (x){neg.u.top[[x]][neg.u.top[[x]]$ipro_des=="-",]})
#pos.u.blank <- lapply(1:length(pos.u.top), function(x){neg.u.top[[x]][neg.u.top[[x]]$ipro_des=="-",]})
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

##It was stopping every time there was an empty list so I made it to pick datasets that included some entries
#lapply(c(1:length(neg.u.top))[sapply(1:length(neg.u.top),function(y){nrow(neg.u.top[[y]])>0})],function(x){
#    jpeg(paste0(names(type.re[x+1]), "_top30_proteins_unmatched.jpg"), width = 900, height = 900)
#    par(mar=c(4,20,4,4))
#    barplot(head(sort(table(neg.u.top[[x]]$ipro_des), decreasing=T),30), horiz=T, type="l",  main=paste0(names(type.re[x+1])\
," top 30 counts"), las=2, cex.names=0.75)
#    box()
#    dev.off()
#})

###A venn diagram function that can compare between two different datasets in a list ( in this case each element of the list \
are the matches in the different dataset)
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

#For pfam and gene3d
venn.d(1,5,type="accession", data.l=neg.u.data)
#For Pfam and superfamily
venn.d(1,6)

###Some specific numbers

lapply(c(6,5,2,1),
       function(x){
           par(mar=c(5,8,8,2))
           common <- intersect(neg.u.top[[x]]$ipro_des, neg.m.top[[x]]$ipro_des)
           un.uniq <- subset(neg.u.top[[x]], !(neg.u.top[[x]]$ipro_des%in% common))
           head(sort(table(un.uniq$ipro_des), decreasing=T),10)
       })
###Or a better look for the unmatching/matching ones
p.db <- 6
common <- intersect(neg.u.top[[p.db]]$ipro_des, neg.m.top[[p.db]]$ipro_des)
c("matched_in_common",nrow(subset(neg.m.top[[p.db]], (neg.m.top[[p.db]]$ipro_des %in% common))))
c("matched_unique",nrow(subset(neg.m.top[[p.db]], !(neg.m.top[[p.db]]$ipro_des %in% common))))
c("unmatched_in_common",nrow(subset(neg.u.top[[p.db]], (neg.u.top[[p.db]]$ipro_des %in% common))))
c("unmatched_unique",nrow(subset(neg.u.top[[p.db]], !(neg.u.top[[p.db]]$ipro_des %in% common))))

##Making a scatterplot out of unmatching and matching values
##saved in functions file
col.a <- colorRampPalette(c( "skyblue", "navy", "violet","orangered"))
col.a <- colorRampPalette(c( "navy", "dodgerblue", "purple","hotpink","red"))

scatplot.proteins(1)

###
###
###make scatterplots using ggplot
##saved in functions file

ggplot.prot.acc.des(1)
dev.off()

ggplot.prot <- function(x, gtype="ipro_des"){
    grid.arrange(scatplot.proteins.gg(x, neg.m.data, neg.u.data, gtype=gtype, mt=paste0(names(type.re[x+1]), ", all transcrib\
ed queries")),
                 scatplot.proteins.gg(x, gtype=gtype, mt=paste0(names(type.re[x+1]),", high scoring matches (evalue-10)")),
                 ##                 scatplot.proteins.gg(x, neg.m.top.nonull, neg.u.top.nonull, gtype=gtype, mt=paste0(names(\
type.re[x+1]),", high scoring matches (evalue-10), exclude unidentified matches")),
                 ncol=2)
}

sapply(1:length(db.bit), function(x){
    jpeg(paste0(names(db.bit[x]),"_protein_scatpl_ipro_des.jpg"), width = 900, height = 1260)
    ggplot.prot(x)
    dev.off()
    })

ggplot.prot.acc.des <- function(x){
    grid.arrange(scatplot.prot.acc.des(x, neg.m.data, neg.u.data, mt=paste0(names(type.re[x+1]), ", all transcribed queries")\
),
                 scatplot.prot.acc.des(x, mt=paste0(names(type.re[x+1]),", high scoring matches (evalue-10)")),

                 ncol=2)
}

scatplot.prot.agg.des(1)
ggplot.prot(2, gtype="description")
ggplot.prot(2, gtype="signature")

