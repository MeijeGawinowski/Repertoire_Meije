# gamma
file1="data-simul-chain3/D1001DT0.06T.fasta"
PMST=NULL
ID=str_split(file1,"D")[[1]][2]
time=str_split(file1,"T")[[1]][2]
## example for one fasta file
subdna<-read.dna(file1,format="fasta")
uniqseq <- dna2uniqSequences(subdna)
## get the counts for each sequence
IDcounts<-do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
IDcounts<-as.data.frame(IDcounts[order(-IDcounts[,1]),])
colnames(IDcounts) <- c( 'count')
seqindex<-match(rownames(IDcounts), labels(uniqseq@uniqdna))
## reorder the DNAbin accordingly
ordereddna<-uniqseq@uniqdna[seqindex, ]
uniqdist<-dist.dna(ordereddna,model="raw", as.matrix=TRUE)
mstdist<-mst(uniqdist)
plotcord <- data.frame(gplot.layout.fruchtermanreingold(mstdist, NULL))
X1=X2=Y1=Y2=NULL
colnames(plotcord) = c("X1","X2")
rownames(plotcord) = rownames(uniqdist)
plotcord<-cbind(plotcord,IDcounts)
mstdist[lower.tri(mstdist,diag=TRUE)]<-NA
eList <- NULL
for ( i in 1:nrow(mstdist) ){
  for ( j in 1:ncol(mstdist)) {
    if (!is.na(mstdist[i,j])){
      if (mstdist[i,j]>0){
        eList <- rbind(eList,c( rownames(mstdist)[i], colnames(mstdist)[j]))
      }
    }
  }
}
edges <- data.frame(plotcord[eList[,1],1:2], plotcord[eList[,2],1:2])
colnames(edges) <-  c("X1","Y1","X2","Y2")
myTheme <- theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),	panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(), panel.border = element_blank(),
                 panel.background = element_blank(), legend.position = "none")
pmst1 <- ggplot()+theme_void()
pmst1 <- pmst1 + geom_segment(data=edges,aes(x=X1,xend=X2,y=Y1,yend=Y2),lineend="round",size=0.5)
pmst1 <- pmst1 + scale_y_continuous("",labels=NULL)+scale_x_continuous("",labels=NULL)
pmst1 <- pmst1 + geom_point(aes(X1, X2, size=count, colour="red"), data=plotcord,colour="turquoise4")
pmst1 <- pmst1 + ggtitle(paste("ID",ID,", t=",time,", With coefficients"))
pmst1




# gamma
file2="data-simul-chain3/D1001DT8.88T.fasta"
PMST=NULL
ID=str_split(file2,"D")[[1]][2]
time=str_split(file2,"T")[[1]][2]
## example for one fasta file
subdna<-read.dna(file2,format="fasta")
uniqseq <- dna2uniqSequences(subdna)
## get the counts for each sequence
IDcounts<-do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
IDcounts<-as.data.frame(IDcounts[order(-IDcounts[,1]),])
colnames(IDcounts) <- c( 'count')
seqindex<-match(rownames(IDcounts), labels(uniqseq@uniqdna))
## reorder the DNAbin accordingly
ordereddna<-uniqseq@uniqdna[seqindex, ]
uniqdist<-dist.dna(ordereddna,model="raw", as.matrix=TRUE)
mstdist<-mst(uniqdist)
plotcord <- data.frame(gplot.layout.fruchtermanreingold(mstdist, NULL))
X1=X2=Y1=Y2=NULL
colnames(plotcord) = c("X1","X2")
rownames(plotcord) = rownames(uniqdist)
plotcord<-cbind(plotcord,IDcounts)
mstdist[lower.tri(mstdist,diag=TRUE)]<-NA
eList <- NULL
for ( i in 1:nrow(mstdist) ){
  for ( j in 1:ncol(mstdist)) {
    if (!is.na(mstdist[i,j])){
      if (mstdist[i,j]>0){
        eList <- rbind(eList,c( rownames(mstdist)[i], colnames(mstdist)[j]))
      }
    }
  }
}
edges <- data.frame(plotcord[eList[,1],1:2], plotcord[eList[,2],1:2])
colnames(edges) <-  c("X1","Y1","X2","Y2")
myTheme <- theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),	panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(), panel.border = element_blank(),
                 panel.background = element_blank(), legend.position = "none")
pmst2 <- ggplot()+theme_void()
pmst2 <- pmst2 + geom_segment(data=edges,aes(x=X1,xend=X2,y=Y1,yend=Y2),lineend="round",size=0.45)
pmst2 <- pmst2 + scale_y_continuous("",labels=NULL)+scale_x_continuous("",labels=NULL)
pmst2 <- pmst2 + geom_point(aes(X1, X2, size=count, colour="red"), data=plotcord,colour="turquoise4")
pmst2 <- pmst2 + ggtitle(paste("ID",ID,", t=",time,", With coefficients"))
pmst2

# no gamma
file3="data-simul-chain4/D1001DT0.57T.fasta"
PMST=NULL
ID=str_split(file3,"D")[[1]][2]
time=str_split(file3,"T")[[1]][2]
## example for one fasta file
subdna<-read.dna(file3,format="fasta")
uniqseq<- dna2uniqSequences(subdna)
## get the counts for each sequence
IDcounts<-do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
IDcounts<-as.data.frame(IDcounts[order(-IDcounts[,1]),])
colnames(IDcounts) <- c( 'count')
seqindex<-match(rownames(IDcounts), labels(uniqseq@uniqdna))
## reorder the DNAbin accordingly
ordereddna<-uniqseq@uniqdna[seqindex, ]
uniqdist<-dist.dna(ordereddna,model="raw", as.matrix=TRUE)
mstdist<-mst(uniqdist)
plotcord <- data.frame(gplot.layout.fruchtermanreingold(mstdist, NULL))
X1=X2=Y1=Y2=NULL
colnames(plotcord) = c("X1","X2")
rownames(plotcord) = rownames(uniqdist)
plotcord<-cbind(plotcord,IDcounts)
mstdist[lower.tri(mstdist,diag=TRUE)]<-NA
eList <- NULL
for ( i in 1:nrow(mstdist) ){
  for ( j in 1:ncol(mstdist)) {
    if (!is.na(mstdist[i,j])){
      if (mstdist[i,j]>0){
        eList <- rbind(eList,c( rownames(mstdist)[i], colnames(mstdist)[j]))
      }
    }
  }
}
edges <- data.frame(plotcord[eList[,1],1:2], plotcord[eList[,2],1:2])
colnames(edges) <-  c("X1","Y1","X2","Y2")
myTheme <- theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),	panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(), panel.border = element_blank(),
                 panel.background = element_blank(), legend.position = "none")
pmst3 <- ggplot()+theme_void()
pmst3 <- pmst3 + geom_segment(data=edges,aes(x=X1,xend=X2,y=Y1,yend=Y2),lineend="round",size=0.45)
pmst3 <- pmst3 + scale_y_continuous("",labels=NULL)+scale_x_continuous("",labels=NULL)
pmst3 <- pmst3 + geom_point(aes(X1, X2, size=count, colour="red"), data=plotcord,colour="turquoise4")
pmst3 <- pmst3 + ggtitle(paste("ID",ID,", t=",time,", Without coefficients"))
pmst3

file4="data-simul-chain2/D1001DT8.88T.fasta"
PMST=NULL
ID=str_split(file4,"D")[[1]][2]
time=str_split(file4,"T")[[1]][2]
## example for one fasta file
subdna<-read.dna(file1,format="fasta")
uniqseq <- dna2uniqSequences(subdna)
## get the counts for each sequence
IDcounts<-do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
IDcounts<-as.data.frame(IDcounts[order(-IDcounts[,1]),])
colnames(IDcounts) <- c( 'count')
seqindex<-match(rownames(IDcounts), labels(uniqseq@uniqdna))
## reorder the DNAbin accordingly
ordereddna<-uniqseq@uniqdna[seqindex, ]
uniqdist<-dist.dna(ordereddna,model="raw", as.matrix=TRUE)
mstdist<-mst(uniqdist)
plotcord <- data.frame(gplot.layout.fruchtermanreingold(mstdist, NULL))
X1=X2=Y1=Y2=NULL
colnames(plotcord) = c("X1","X2")
rownames(plotcord) = rownames(uniqdist)
plotcord<-cbind(plotcord,IDcounts)
mstdist[lower.tri(mstdist,diag=TRUE)]<-NA
eList <- NULL
for ( i in 1:nrow(mstdist) ){
  for ( j in 1:ncol(mstdist)) {
    if (!is.na(mstdist[i,j])){
      if (mstdist[i,j]>0){
        eList <- rbind(eList,c( rownames(mstdist)[i], colnames(mstdist)[j]))
      }
    }
  }
}
edges <- data.frame(plotcord[eList[,1],1:2], plotcord[eList[,2],1:2])
colnames(edges) <-  c("X1","Y1","X2","Y2")
myTheme <- theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),	panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(), panel.border = element_blank(),
                 panel.background = element_blank(), legend.position = "none")
pmst4 <- ggplot()+theme_void()
pmst4 <- pmst4 + geom_segment(data=edges,aes(x=X1,xend=X2,y=Y1,yend=Y2),lineend="round",size=0.45)
pmst4 <- pmst4 + scale_y_continuous("",labels=NULL)+scale_x_continuous("",labels=NULL)
pmst4 <- pmst4 + geom_point(aes(X1, X2, size=count, colour="red"), data=plotcord,colour="turquoise4")
pmst4 <- pmst4 + ggtitle(paste("ID",ID,", t=",time,", Without coefficients"))
pmst4

PMST=c(list(pmst1),list(pmst2))
do.call(grid.arrange,PMST)
PMST2=c(list(pmst1),list(pmst2),list(pmst3),list(pmst4))
do.call(grid.arrange,PMST2)
PMST3=c(list(pmst1),list(pmst3))
do.call(PMST3)
