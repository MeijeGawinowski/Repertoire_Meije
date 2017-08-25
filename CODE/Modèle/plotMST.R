if(FALSE){

    library(ape)
    library(ggplot2)
    library(sna)
    library(gridExtra)
    library(stringr)
    library(grDevices)

}

## Function to convert a class 'DNAbin' to the class 'uniqSequences'
## param dna an object of the class "DNAbin"
## author Joseph Hughes
dna2uniqSequences<-function(dna){
    if(is.null(dna)) return(NULL)
    if(!inherits(dna,"DNAbin")) {
        warning("dna should be a DNAbin object - returning NULL")
        return(NULL)
    }
    seqmatrix <- as.character(dna)
    seqstring<-apply(format(seqmatrix), 1, paste, collapse="")
    uniqList<-tapply(names(seqstring),list(seqstring),I)
    splitseqstr<-strsplit(names(uniqList),"")
    uniqmat<-do.call(rbind, splitseqstr)
    uniqnames<-paste("uniqseqID",seq(1:length(uniqList)),sep="")
    rownames(uniqmat)<-uniqnames
    names(uniqList)<-uniqnames
    uniqdna<-as.DNAbin(uniqmat)
    uniqSequences<-new("uniqSequences",uniqID=uniqList,uniqdna=uniqdna)


  return(uniqSequences)
}


############################
####  CLASSE DEFINITION ####
############################

## CLASS DESCRIPTION:
## Instance of uniqSequences store sequences; its content includes:
## - @uniqID: data about the identical sequences, stored as a list of vectors of sequenceID
## - @uniqdna: unique dna sequences, stored as a DNAbin

## setOldClass("DNAbin")
setClass("uniqSequences", representation(uniqID="list", uniqdna="DNAbin"),prototype(uniqID=NULL, uniqdna=NULL))


######################
####  CONSTRUCTOR ####
######################

## INPUT DESCRIPTION:
## 'uniqID': a list of vectors, each vector contains the sequenceID of sequences that are identical
## 'uniqdna': a DNAbin with named unique sequences according to the list names
##
setMethod("initialize", "uniqSequences", function(.Object, uniqID=NULL, uniqdna=NULL){
    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object

    ## store old option ##
    o.opt <- options("stringsAsFactors")
    options("stringsAsFactors"=FALSE)
    on.exit(options(o.opt))

    ## PROCESS INFORMATION ABOUT uniqIDs ('uniqID') ##
    if (!is.null(uniqdna) && inherits(uniqdna,"DNAbin")){
      x@uniqdna<-uniqdna
    }
    if(!is.null(uniqID)){
      x@uniqID <- list()
      for(i in 1:length(uniqID))
        {
          x@uniqID[[i]] <- uniqID[[i]]
        }
      names(x@uniqID) <- names(uniqID)
      ## make sure that all the uniqIDs are in 'uniqdna'
      if(!is.null(x@uniqID) && !is.null(x@uniqdna)){
        # print(labels(x@uniqdna))
        # print(names(x@uniqID))
        unknownIDs <- unique(labels(x@uniqdna))[!unique(labels(x@uniqdna)) %in% names(x@uniqID)]
        if(length(unknownIDs)>0) {
          unknownIDs.txt <- paste(unknownIDs, collapse = ", ")
          warning(paste("the following uniqIDs in the DNAbin do not have information in the list:\n", unknownIDs.txt))
        }
        unknownlabs <- unique(labels(x@uniqID))[!unique(names(x@uniqID)) %in% labels(x@uniqdna)]
        if(length(unknownlabs)>0) {
          unknownlabs.txt <- paste(unknownlabs, collapse = ", ")
          warning(paste("the following uniqIDs in the list do not have sequences in the DNAbin:\n", unknownlabs.txt))
        }

      }
    }

    ## RETURN OBJECT ##
    return(x)
}) # end uniqSequences constructor



###############################################################

filenames=paste(directory,system(paste("ls",directory),intern=TRUE),
    sep="")
filenames=filenames[-length(filenames)]
PMST=NULL

#subdna=read.dna(filenames[1],format="fasta")
#for(ff in filenames[-1]){
## example for one fasta file
#subdna<-rbind(subdna,read.dna(ff,format="fasta"))
#}

for(ff in filenames){
    obsID <- str_split(ff,"D")[[1]][2]
    obstime <- str_split(ff,"T")[[1]][2]
    ## example for one fasta file
    subdna<-read.dna(ff,format="fasta")
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
    pmst <- ggplot() + theme_classic()
    pmst <- pmst + geom_segment(data=edges,aes(x=X1,xend=X2,y=Y1,yend=Y2),lineend="round",size=0.5)
    pmst <- pmst + scale_y_continuous("",labels=NULL)+scale_x_continuous("",labels=NULL)
    pmst <- pmst + geom_point(aes(X1, X2, size=count, colour="Blue"), data=plotcord,colour="turquoise4")
    pmst <- pmst + ggtitle(paste("ID ",obsID,", t = ",obstime))
    pmst
    PMST=c(PMST,list(pmst))
}


do.call(grid.arrange,PMST)



