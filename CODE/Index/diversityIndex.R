rm(list=ls())
library(ape)
library(vegan)
library(pegas)
library(poppr)
library(DiversitySampler)
library(untb)
library(reshape2)
library(ggplot2)

# function that counts the number of non null elements in a vector
countProp=function(v){
  c=0
  for (i in 1:length(v)){
    if (v[i] != 0){c=c+1}
  }
  return(c)
}


# function that turns a sample into DNAbin sequences
toDNAbin=function(samp,sequence){
  # samp : indexes of the sampled sequences
  # sequence : matrix with all the sequences
  M=dim(sequence)[2]
  v=matrix(0,nrow=length(samp),ncol=M)
  for (i in 1:length(samp)){
    v[i,]=sequence[samp[i],]
  }
  v=c(v)
  v=as.character(v)
  for (i in 1:length(v)){
    if(v[i]=="1"){v[i]="A"}
    if(v[i]=="2"){v[i]="C"}
    if(v[i]=="3"){v[i]="G"}
    if(v[i]=="4"){v[i]="T"}
  }
  v=as.DNAbin(matrix(c(v),ncol=M,nrow=length(samp)))
  return(v)
}

## function that describes the viral kinetics in one host
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


# function that assesses diversity at each moment, returns a data frame with time and index
intraHostDiversity=function(init=list(prop=1,sequ1=rep(1,10^4)),kinetics,mu,gamma,fitnessRule,method=1,size=2){
  ## mu: mutation rate per nucleotide, per virus particle (i.e. virion), per unit of time
  ## gamma : evolution rates (effects of the environment)
  ## tolerance: discard variants with proportions less than tolerance
  ## fitnessRule: function determining how fit is a variant
  if(nrow(kinetics)==1){browser()}
  M=length(init[[2]])
  lags=diff(kinetics$time)
  nmutation=0
  ## initialize the composition of variants
  compo=init
  ## number of information elements provided in compo before the sequences (here, only 1 : prop)
  ninfo=1
  ## compute fitness for initial variants
  fitness=NULL
  
  for(k in (ninfo+1):length(compo)){
    fitness=c(fitness,fitnessRule(compo[[k]]))
  }
  
  nb_vector=c()
  Richness=c()
  Time=c()
  Shannon=c()
  indice=c()
  
  for(i in 2:length(kinetics$time)){
    if (kinetics$virus[i] >= 1){
      print(i)
      Time=c(Time,i)      
      ## draw the number of mutations for time step i
      nmut=rpois(1,mu*M*kinetics$virus[i-1]*lags[i-1])
      nmutation=c(nmutation,nmut)
      if(nmut>0){
        propCand=rep(1,nmut)
        ## draw which variants are used for generating new variants
        mutatingSequ=sample(x=1:(length(compo)-ninfo),size=length(propCand),
                            replace=TRUE,prob=compo$prop)
        ## update the vector of proportions
        compo$prop=c(compo$prop*(1-nmut/kinetics$virus[i-1]),
                     propCand*nmut/kinetics$virus[i-1])
        
        ## draw the mutations
        for(k in 1:length(propCand)){
          newsequ=compo[[ninfo+mutatingSequ[k]]]
          mutatingNucl=sample(1:M,1)
          newsequ[mutatingNucl]=
            sample((1:4)[(1:4)!=newsequ[mutatingNucl]],1)
          compo=c(compo,list(newsequ))
          ## compute fitness of new variants
          fitness=c(fitness,fitnessRule(newsequ))
        }
        
      }
      ## multiplication probability depends on fitness
      pmult=compo$prop*fitness
      pmult=pmult/sum(pmult)
      pmult=rnorm(length(pmult),mean=pmult,sd=sqrt(gamma[1]*pmult*(1-pmult)^gamma[2]))
      pmult=pmin(1,pmax(0,pmult))
      ## check that pmult is a probability
  
      compo$prop=as.numeric(rmultinom(n=1,size=kinetics$virus[i],
                                      prob=pmult))

      compo$prop=compo$prop/sum(compo$prop)

      ## compute new proportions of variants
      proportion=compo$prop
      N=length(proportion)
      sequence=matrix(0,N,M)
      for(k in 1:N){
        sequence[k,]=compo[[ninfo+k]]
      }
      
      ## Richness estimator
      if (method==1){
        ylab = "Richness estimator"
        title="data_richness.csv"
        indice=c(indice,countProp(compo$prop))
      }
      
      ## Shannon Diversity Index
      if (method==2){
        ylab="Shannon Index"
        title="data_shannon.csv"
        nb=compo$prop*kinetics$virus[i]
        H1=diversity(nb,index="shannon",base=exp(1))
        indice=c(indice,H1)
      }
      
      ## Simpson Index
      if (method==3){
        ylab="Simpson Index"
        title="data_simpson.csv"
        nb=compo$prop*kinetics$virus[i]
        H2=diversity(nb,index="simpson")
        indice=c(indice,H2)
      }
      
      ## Pielou's Evenness
      if(method==4){
        ylab="Pielou Evenness"
        title="data_pielou.csv"
        Richness=c(Richness,countProp(compo$prop))
        nb=compo$prop*kinetics$virus[i]
        Shannon=c(Shannon,diversity(nb,index="shannon",base=exp(1)))
        indice=Shannon/log10(Richness)
      }
      
      ## Nucleotide Diversity or p-distance
      if(method==5){
        ylab="p-distance"
        title="data_nuc.csv"
        samp=sample(x=1:(length(compo)-ninfo),size=size,replace=TRUE,prob=compo$prop)
        seq=toDNAbin(samp,sequence)
        indice=c(indice,nuc.div(seq))
      }
      
      ## Jukes and Cantor Distance
      if(method==6){
        ylab="Jukes and Cantor distance"
        title="data_JC.csv"
        samp=sample(x=1:(length(compo)-ninfo),size=size,replace=TRUE,prob=compo$prop)
        seq=toDNAbin(samp,sequence)
        indice=c(indice,mean(dist.dna(seq,model="JC69")))
      }
      
      ## Nei Distance
      if(method==7){
        ylab="Nei distance"
        title="data_nei.csv"
        samp=sample(x=1:(length(compo)-ninfo),size=size,replace=TRUE,prob=compo$prop)
        seq=toDNAbin(samp,sequence)
        indice=c(indice,mean(nei.dist(as.genind(seq))))
      }
      
    }
  }
  data=data.frame(Time,indice)
  write.table(data,title,sep="\t")
  return(list(data=data,ylab=ylab))
}


## PARAMETERS
vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=5.0e-05,delta=2.0,p=1e-06,c=0.5,timelag=0.1,duration=45)
kin1 = viralKinetics(init=vpar1$init,vpar1$s,vpar1$d,vpar1$beta,vpar1$delta,vpar1$p,vpar1$c,vpar1$timelag,vpar1$duration)
spar1=list(init=list(prop=1,sequ1=rep(1,1000)),mu=5e-06,gamma=c(0.001,1000),tolerance=1e-05)
spar2=list(init=list(prop=1,sequ1=rep(1,1000)),mu=5e-06,gamma=c(0,1000),tolerance=1e-05)
fpar2=function(sequ){
  M=length(sequ)
  return(1*(1+9*mean(sequ[1:round(0.1*M)]==2))/
           (1+9*mean(sequ[1:round(0.1*M)]==3)))
}
fpar1=function(sequ){ return(1) }

## gamma, no fitness rule -> spar1, fpar1 (REFERENCE)
setwd("~/Ref_files")
richness1=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,method=1)
shannon1=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,method=2)
simpson1=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,method=3)
pielou1=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,method=4)
nucdiv1=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,method=5,size=100)
JCdistance1=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,method=6,size=100)
neidistance1=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,method=7,size=100)

## no gamma, no fitness rule -> spar2, fpar1
setwd("~/gamma_files")
richness2=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar2$gamma,fitnessRule=fpar1,method=1)
shannon2=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar2$gamma,fitnessRule=fpar1,method=2)
simpson2=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar2$gamma,fitnessRule=fpar1,method=3)
pielou2=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar2$gamma,fitnessRule=fpar1,method=4)
nucdiv2=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar2$gamma,fitnessRule=fpar1,method=5)
JCdistance2=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar2$gamma,fitnessRule=fpar1,method=6,size=100)
neidistance2=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar2$gamma,fitnessRule=fpar1,method=7,size=100)

## gamma, fitness rule -> spar1, fpar2
setwd("~/fitness_files")
richness3=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar2,method=1)
shannon3=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar2,method=2)
simpson3=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar2,method=3)
pielou3=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar2,method=4)
nucdiv3=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar2,method=5,size=100)
JCdistance3=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar2,method=6,size=100)
neidistance3=intraHostDiversity(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar2,method=7,size=100)







