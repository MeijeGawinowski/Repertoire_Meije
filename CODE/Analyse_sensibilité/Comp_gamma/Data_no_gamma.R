rm(list=ls())
library(ape)
library(vegan)
library(pegas)
library(poppr)
library(DiversitySampler)
library(untb)


## function that turns a sample into DNAbin sequences
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

## fitness rule function
fpar1=function(sequ){ 
  return(1) 
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


# function that assesses diversity at each moment
intraHostDiversity2=function(init=list(prop=1,sequ1=rep(1,10^4)),kinetics,mu,gamma,fitnessRule,size=2){
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
  
  for(i in 2:length(kinetics$time)){
    if (kinetics$virus[i] >= 1){
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
      
    }
  }
  ## Nucleotide Diversity or p-distance
  samp=sample(x=1:(length(compo)-ninfo),size=size,replace=TRUE,prob=compo$prop)
  seq=toDNAbin(samp,sequence)
  indice=nuc.div(seq)
  return(indice)
}

## function that runs sensitivity analysis
runSensitivity=function(tf,par,min,max,lag){
  i=min # min value of the tested parameter
  v_time=c()
  v_index=c()
  while (i <= max){
    print(i)
    v_time=c(v_time,i)
    if (par=="beta"){
      vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=i,delta=2.0,p=1e-06,c=0.5,timelag=0.1,duration=tf)
      spar1=list(init=list(prop=1,sequ1=rep(1,1000)),mu=5e-06,gamma=c(0,1000),tolerance=1e-05) 
    }
    if (par=="delta"){
      vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=5.0e-05,delta=i,p=1e-06,c=0.5,timelag=0.1,duration=tf)
      spar1=list(init=list(prop=1,sequ1=rep(1,1000)),mu=5e-06,gamma=c(0,1000),tolerance=1e-05) 
    }
    if (par=="p"){
      vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=5.0e-05,delta=2.0,p=i,c=0.5,timelag=0.1,duration=tf)
      spar1=list(init=list(prop=1,sequ1=rep(1,1000)),mu=5e-06,gamma=c(0,1000),tolerance=1e-05) 
    }
    if (par=="c"){
      vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=5.0e-05,delta=2.0,p=1e-06,c=i,timelag=0.1,duration=tf)
      spar1=list(init=list(prop=1,sequ1=rep(1,1000)),mu=5e-06,gamma=c(0,1000),tolerance=1e-05) 
    }
    if (par=="mu"){
      vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=5.0e-05,delta=2.0,p=1e-06,c=0.5,timelag=0.1,duration=tf)
      spar1=list(init=list(prop=1,sequ1=rep(1,1000)),mu=i,gamma=c(0,1000),tolerance=1e-05) 
    }
    
    kin1 = viralKinetics(init=vpar1$init,vpar1$s,vpar1$d,vpar1$beta,vpar1$delta,vpar1$p,vpar1$c,vpar1$timelag,vpar1$duration)
    v_index=c(v_index,intraHostDiversity2(init=spar1$init,kinetics=kin1,mu=spar1$mu,gamma=spar1$gamma,fitnessRule=fpar1,size=100))
    
    i=i+lag
  }
  return(list(t=v_time,indice=v_index))
}

## WITHIN-HOST DIVERSITY ASSESSMENT (NO FITNESS RULE)
spar2=list(init=list(prop=1,sequ1=rep(1,1000)),mu=5e-06,gamma=c(0,1000),tolerance=1e-05)
fpar2=function(sequ){
  M=length(sequ)
  return(1*(1+9*mean(sequ[1:round(0.1*M)]==2))/
           (1+9*mean(sequ[1:round(0.1*M)]==3)))
}

setwd("~/Comp_gamma/data_gamma")

## beta test for t=3
beta1=runSensitivity(tf=3,par="beta",min=1e-06,max=1e-04,lag=1e-06)
beta2=runSensitivity(tf=3,par="beta",min=1e-06,max=1e-04,lag=1e-06)
beta3=runSensitivity(tf=3,par="beta",min=1e-06,max=1e-04,lag=1e-06)
beta4=runSensitivity(tf=3,par="beta",min=1e-06,max=1e-04,lag=1e-06)
beta5=runSensitivity(tf=3,par="beta",min=1e-06,max=1e-04,lag=1e-06)
beta_t=c(beta1$t,beta2$t,beta3$t,beta4$t,beta5$t)
beta_i=c(beta1$indice,beta2$indice,beta3$indice,beta4$indice,beta5$indice)
data_beta=data.frame(beta_t,beta_i)
colnames(data_beta)=c("par","dist")
write.table(data_beta,"data_beta.csv",sep=",")

## delta test for t=3
delta1=runSensitivity(tf=3,par="delta",min=1,max=14,lag=0.1)
delta2=runSensitivity(tf=3,par="delta",min=1,max=14,lag=0.1)
delta3=runSensitivity(tf=3,par="delta",min=1,max=14,lag=0.1)
delta4=runSensitivity(tf=3,par="delta",min=1,max=14,lag=0.1)
delta5=runSensitivity(tf=3,par="delta",min=1,max=14,lag=0.1)
delta_t=c(delta1$t,delta2$t,delta3$t,delta4$t,delta5$t)
delta_i=c(delta1$indice,delta2$indice,delta3$indice,delta4$indice,delta5$indice)
data_delta=data.frame(delta_t,delta_i)
colnames(data_delta)=c("par","dist")
write.table(data_delta,"data_delta.csv",sep=",")

## p test for t=3
p1=runSensitivity(tf=3,par="p",min=1e-08,max=5e-07,lag=1e-08)
p2=runSensitivity(tf=3,par="p",min=1e-08,max=5e-07,lag=1e-08)
p3=runSensitivity(tf=3,par="p",min=1e-08,max=5e-07,lag=1e-08)
p4=runSensitivity(tf=3,par="p",min=1e-08,max=5e-07,lag=1e-08)
p5=runSensitivity(tf=3,par="p",min=1e-08,max=5e-07,lag=1e-08)
p_t=c(p1$t,p2$t,p3$t,p4$t,p5$t)
p_i=c(p1$indice,p2$indice,p3$indice,p4$indice,p5$indice)
data_p=data.frame(p_t,p_i)
colnames(data_p)=c("par","dist")
write.table(data_p,"data_p.csv",sep=",")

## c test for t=3
c1=runSensitivity(tf=3,par="c",min=0.1,max=1,lag=0.01)
c2=runSensitivity(tf=3,par="c",min=0.1,max=1,lag=0.01)
c3=runSensitivity(tf=3,par="c",min=0.1,max=1,lag=0.01)
c4=runSensitivity(tf=3,par="c",min=0.1,max=1,lag=0.01)
c5=runSensitivity(tf=3,par="c",min=0.1,max=1,lag=0.01)
c_t=c(c1$t,c2$t,c3$t,c4$t,c5$t)
c_i=c(c1$indice,c2$indice,c3$indice,c4$indice,c5$indice)
data_c=data.frame(c_t,c_i)
colnames(data_c)=c("par","dist")
write.table(data_c,"data_c.csv",sep=",")

## mu test for t=3
m1=runSensitivity(tf=3,par="mu",min=1e-07,max=1e-05,lag=1e-07)
m2=runSensitivity(tf=3,par="mu",min=1e-07,max=1e-05,lag=1e-07)
m3=runSensitivity(tf=3,par="mu",min=1e-07,max=1e-05,lag=1e-07)
m4=runSensitivity(tf=3,par="mu",min=1e-07,max=1e-05,lag=1e-07)
m5=runSensitivity(tf=3,par="mu",min=1e-07,max=1e-05,lag=1e-07)
m_t=c(m1$t,m2$t,m3$t,m4$t,m5$t)
m_i=c(m1$indice,m2$indice,m3$indice,m4$indice,m5$indice)
data_m=data.frame(m_t,m_i)
colnames(data_m)=c("par","dist")
write.table(data_m,"data_m.csv",sep=",")
