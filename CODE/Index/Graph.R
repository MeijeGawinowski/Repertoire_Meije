rm(list=ls())

library(reshape2)
library(ggplot2)

scatter_plot=function(data,ylab,colour){
  plot(data$Time,data$indice,ylab=ylab,xlab="time",col=colour)
  return(data)
}

smooth_plot=function(data1,data2,comp,index){
  # file1 must be the reference (gamma, no fitnessı)
  newdata=data.frame(data1$Time,data1$indice,data2$indice)
  colnames(newdata)
  
  if (comp=="Gamma"){
    colnames(newdata)=c("time","Gamma","No Gamma")
  }
  if (comp=="fitness"){
    colnames(newdata)=c("time","No fitness","Fitness")
  }
  newdata=melt(newdata,id.vars="time")
  pl=ggplot(newdata,aes(x=time,y=value,colour=variable))+geom_smooth(aes(fill=variable),method="loess")+theme_classic()+ylab(index)
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


## Files loading
setwd("~/Ref_files")
richness1=read.table("data_richness.csv",sep="\t",header=TRUE)
shannon1=read.table("data_shannon.csv",sep="\t",header=TRUE)
simpson1=read.table("data_simpson.csv",sep="\t",header=TRUE)
pielou1=read.table("data_pielou.csv",sep="\t",header=TRUE)
nuc1=read.table("data_nuc.csv",sep="\t",header=TRUE)
jc1=read.table("data_JC.csv",sep="\t",header=TRUE)
nei1=read.table("data_nei.csv",sep="\t",header=TRUE)
setwd("~/gamma_files")
richness2=read.table("data_richness.csv",sep="\t",header=TRUE)
shannon2=read.table("data_shannon.csv",sep="\t",header=TRUE)
simpson2=read.table("data_simpson.csv",sep="\t",header=TRUE)
pielou2=read.table("data_pielou.csv",sep="\t",header=TRUE)
nuc2=read.table("data_nuc.csv",sep="\t",header=TRUE)
jc2=read.table("data_JC.csv",sep="\t",header=TRUE)
nei2=read.table("data_nei.csv",sep="\t",header=TRUE)
setwd("~/fitness_files")
richness3=read.table("data_richness.csv",sep="\t",header=TRUE)
shannon3=read.table("data_shannon.csv",sep="\t",header=TRUE)
simpson3=read.table("data_simpson.csv",sep="\t",header=TRUE)
pielou3=read.table("data_pielou.csv",sep="\t",header=TRUE)
nuc3=read.table("data_nuc.csv",sep="\t",header=TRUE)
jc3=read.table("data_JC.csv",sep="\t",header=TRUE)
nei3=read.table("data_nei.csv",sep="\t",header=TRUE)

## Viral kinetics
vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=5.0e-05,delta=2.0,p=1e-06,c=0.5,timelag=0.01,duration=50)
viralKinetics=function(init,s,d,beta,delta,p,c,timelag,duration){
  ## MODEL COMPONENTS
  ## DOI: 10.1002/wsbm.129
  ## http://onlinelibrary.wiley.com/doi/10.1002/wsbm.129/full
  ## T: quantity of target cells
  ## I: quantity of infected cells2
  ## V: quantity of infectious viral particles
  ## s: supply rate of target cells
  ## d: death rate of target cells
  ## beta: infection rate
  ## delta: death rate of infected cells
  ## p: production rate of virus
  ## c: clearance rate of virus
  ## OTHER NOTATIONS 
  ## init: initial state for T, I and V
  ## timelag: numerical time lag for solving the equation system
  ## duration: duration of infection
  target=init$TT
  infected=init$I
  virus=init$V
  time=seq(0,duration,by=timelag)
  for(i in 1:(length(time)-1)){
    ## implicit method for ode resolution
    target=c(target,target[i]/(1-timelag*(s-d)+timelag*beta*virus[i]))
    infected=c(infected,(timelag*beta*target[i]*virus[i]+infected[i])/
                 (1+timelag*delta))
    virus=c(virus,(timelag*p*infected[i]+virus[i])/(1+timelag*c))
  }
  ## output: time, nb of target cells, nb of infected cells, nb of virions
  return(data.frame(time,target,infected,virus))
}
kin1=viralKinetics(vpar1$init,vpar1$s,vpar1$d,vpar1$beta,vpar1$delta,vpar1$p,vpar1$c,vpar1$timelag,vpar1$duration)

## Scatter plot
par(mfrow=c(4,2))
plot(kin1$time,kin1$virus,ylab="Virion population",xlab="time")
scatter_plot(richness1,"Richness Estimator",colour="forestgreen")
scatter_plot(shannon1, "Shannon Index",colour="mediumorchid")
scatter_plot(simpson1, "Simpson Index",colour="dodgerblue4")
scatter_plot(pielou1,"Pielou Evenness",colour="deeppink3")
scatter_plot(nuc1,"p-distance",colour="springgreen2")
scatter_plot(jc1,"Jukes and Cantor distance",colour="firebrick2")
scatter_plot(nei1,"Nei distance",colour="steelblue1")

## Smooth plot gamma
pg1=smooth_plot(richness1,richness2,"Gamma",index="Richness Estimator")
pg2=smooth_plot(shannon1,shannon2,"Gamma",index="Shannon Index")
pg3=smooth_plot(simpson1,simpson2,"Gamma",index="Simpson Index")
pg4=smooth_plot(pielou1,pielou2,"Gamma",index="Pielou Evenness")
pg5=smooth_plot(nuc1,nuc2,"Gamma",index="p-distance")
pg6=smooth_plot(jc1,jc2,"Gamma",index="Jukes and Cantor distance")
pg7=smooth_plot(nei1,nei2,"Gamma",index="Nei distance")
multiplot(pg1,pg2,pg3,pg4,pg5,pg6,pg7,cols=2)

## Smooth plot fitness
pf1=smooth_plot(richness1,richness3,"fitness",index="Richness estimator")
pf2=smooth_plot(shannon1,shannon3,"fitness",index="Shannon Index")
pf3=smooth_plot(simpson1,simpson3,"fitness",index="Simpson Index")
pf4=smooth_plot(pielou1,pielou3,"fitness",index="Pielou Evenness")
pf5=smooth_plot(nuc1,nuc3,"fitness",index="p-distance")
pf6=smooth_plot(jc1,jc3,"fitness",index="Jukes and Cantor distance")
pf7=smooth_plot(nei1,nei3,"fitness",index="Nei distance")
multiplot(pf1,pf2,pf3,pf4,pf5,pf6,pf7,cols=2)



