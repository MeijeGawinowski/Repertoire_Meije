rm(list=ls())

## load data from data_ref file
data_beta=read.table("data_beta.csv",sep=",",head=TRUE)
data_delta=read.table("data_delta.csv",sep=",",header=TRUE)
data_p=read.table("data_p.csv",sep=",",header=TRUE)
data_c=read.table("data_c.csv",sep=",",header=TRUE)
data_mu=read.table("data_m.csv",sep=",",header=TRUE)
data_g1=read.table("data_g1.csv",sep=",")
data_g2=read.table("data_g2.csv",sep=",")

## load data from data_fitness file
setwd("~/Comp_fitness/data_fitness")
data_beta_f=read.table("data_beta.csv",sep=",",head=TRUE)
data_delta_f=read.table("data_delta.csv",sep=",",header=TRUE)
data_p_f=read.table("data_p.csv",sep=",",header=TRUE)
data_c_f=read.table("data_c.csv",sep=",",header=TRUE)
data_mu_f=read.table("data_m.csv",sep=",",header=TRUE)
data_g1_f=read.table("data_g1.csv",sep=",")
data_g2_f=read.table("data_g2.csv",sep=",")


Graph_loess=function(data1,data2,main){
  # par=data[,2]
  # div=data[,3]
  # newdata=data.frame(par,div)
  data=data.frame(data1$par,data1$dist,data2$dist)
  colnames(data)=c("Parameter","No fitness","Fitness")
  data=melt(data,id.vars="Parameter")
  pl=ggplot(data=data,aes(x=Parameter,y=value,colour=variable))+geom_smooth(aes(fill=variable),method="loess")+theme_classic() 
  pl=pl+ggtitle(main)+xlab("Parameter")+ylab("p-distance")
  return(pl)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(reshape2)
library(ggplot2)

p_beta=Graph_loess(data_beta,data_beta_f,"Beta parameter")
p_delta=Graph_loess(data_delta,data_delta_f,"Delta parameter")
p_p=Graph_loess(data_p,data_p_f,"p parameter")
p_c=Graph_loess(data_c,data_c_f,"c parameter")
p_m=Graph_loess(data_mu,data_mu_f,"mu parameter")
p_g1=Graph_loess(data_g1,data_g1_f,"gamma1 parameter")
p_g2=Graph_loess(data_g1,data_g1_f,"gamma2 parameter")
multiplot(p_beta,p_delta,p_p,cols=3)
multiplot(p_c,p_m,p_g1,p_g2,cols=2)
multiplot(p_m,p_g1,cols=2)
print(p_g2)
