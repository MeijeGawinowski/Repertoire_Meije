library(ape)
library(ggplot2)
library(sna)
library(gridExtra)

rm(list=ls())

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
        ## implicit method of ode resolution
        target=c(target,target[i]/(1-timelag*(s-d)+timelag*beta*virus[i]))
        infected=c(infected,(timelag*beta*target[i]*virus[i]+infected[i])/
            (1+timelag*delta))
        virus=c(virus,(timelag*p*infected[i]+virus[i])/(1+timelag*c))
    }
    ## output: time, nb of target cells, nb of infected cells, nb of virions
    return(data.frame(time,target,infected,virus))
}


# function that describes the viral genetic composition in one host
viralCompo=function(init=list(prop=1,sequ1=rep(1,10^4)),kinetics,mu,gamma,
    tolerance,fitnessRule){
    ## mu : mutation rate per nucleotide, per virus particle (i.e. virion), per unit of time
    ## gamma : vector with two components of genetic evolution representing the effect of the environment
    ## tolerance : discard variants with proportions less than tolerance
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
        ## addition of a Gaussian noise
        pmult=pmin(1,pmax(0,pmult))
        ## 0-1 inflation on pmult
        compo$prop=as.numeric(rmultinom(n=1,size=kinetics$virus[i],
            prob=pmult))
        compo$prop=compo$prop/sum(compo$prop)
        ## compute new proportions of variants
    }
    
    proportion=compo$prop
    N=length(proportion)
    sequence=matrix(0,N,M)
    
    for(k in 1:N){
        sequence[k,]=compo[[ninfo+k]]
    }
    
    sequence=rbind(sequence[proportion>tolerance,])
    fitness=fitness[proportion>tolerance]
    proportion=proportion[proportion>tolerance]
    proportion=proportion/sum(proportion)
    N=length(proportion)
    alreadyStored=duplicated(sequence)
    
    ## check if some sequences are identicals, it it the case
    ## they're brought together and their proportions are added
    for(k in (1:N)[alreadyStored]){
        firstStored=((1:N)[apply(sequence,1,function(u)
                                 sum(sequence[k,]==u)==M)])[1]
        proportion[firstStored]=proportion[firstStored]+proportion[k]
        proportion[k]=0
    }
    
    sequence=rbind(sequence[proportion>0,])
    fitness=fitness[proportion>0]
    proportion=proportion[proportion>0]
    print(sum(nmutation))
    return(list(composition=list(proportion=proportion,fitness=fitness,
                    sequence=sequence),kinetics=data.frame(kinetics,nmutation)))
    ## output : list with within host viral composition (sequences and their proportions) and the viral kinetics
}



## simulate within-host viral compositions and sample the virus for observation and/or transmission
withinHostData=function(init,kinetics,mu,gamma,tolerance,fitnessRule,observationTime,
    samplesize,directory,identifier,id_chain,transdyn,transmissionRate=NULL){
    ## init : list that contains time and prop
    ## kinetics : data generated by viralKinetics
    ## mu : mutation rate
    ## gamma : evolution rates
    ## tolerance : threshold below which variants are eliminated
    ## fitnessRule : describes how fit a variant is
    ## observationTime : moments chosen by user to observe
    ## samplesize : size of the samples for observation (for the samples of sequences)
    ## identifier : host ID
    ## transmissionTime : moments of transmission (determined by transmit function)
    ## transmissionRate 
  
    oldstate=init
    oldstate$time=NULL
    transmissionTime=transdyn$times
    transmissionHost=transdyn$hosts
    time=unique(sort(c(init$time,observationTime,transmissionTime)))
    transmission=NULL
    
    for(i in 2:length(time)){
        newstate=viralCompo(init=oldstate,
            kinetics=rbind(kinetics[kinetics$time>=time[i-1] &
                kinetics$time<=time[i],]),
            mu=mu,gamma=gamma,tolerance=tolerance,fitnessRule=fitnessRule)
        prop=newstate$composition$proportion
        sequ=newstate$composition$sequence
        
        if(sum(observationTime==time[i])>0){ ## time[i] is an observation time
            ii=(1:length(observationTime))[observationTime==time[i]]
            sampledVariants=sample(x=1:length(prop),size=samplesize[ii],
                replace=TRUE,prob=prop)
            j=0
            
            for(k in sampledVariants){
                j=j+1
                sequACGT=as.character(sequ[k,])
                sequACGT[sequACGT=="1"]="A"
                sequACGT[sequACGT=="2"]="C"
                sequACGT[sequACGT=="3"]="G"
                sequACGT[sequACGT=="4"]="T"
                filename=paste(directory,id_chain,"T",time[i],"TD",identifier,"D.fasta",sep="")
                write.table(paste(">",identifier,"D",time[i],"-variant-",j,
                                  sep=""),file=filename,quote=FALSE,
                            append=(j!=1),col.names=FALSE,row.names=FALSE)
                write.table(paste(sequACGT,collapse=""),file=filename,
                            quote=FALSE,
                            append=TRUE,col.names=FALSE,row.names=FALSE)
            }
            
        } 
        
        else {}
        
        if(sum(transmissionTime==time[i])>0){ ## time[i] is a transmission time
          nexthost=transdyn$hosts[which(transdyn$times==time[i])]
          transmittedVirus=round(kinetics$virus[kinetics$time==time[i]]*
                transmissionRate)
            ## number of transmitted virions to another host
            transmittedProp=as.numeric(rmultinom(n=1,size=transmittedVirus,
                prob=prop))
            transmittedProp=transmittedProp/sum(transmittedProp)
            transmittedSequ=rbind(sequ[transmittedProp>0,])
            transmittedProp=transmittedProp[transmittedProp>0]
            transmittedInit=list(prop=transmittedProp)
            
            for(k in 1:length(transmittedProp)){
                transmittedInit=c(transmittedInit,list(transmittedSequ[k,]))
            }
            
            transmission[[length(transmission)+1]]=
                list(source=identifier,
                     target=nexthost,
                    transmittedTime=time[i],
                     transmittedVirus=transmittedVirus,
                     transmittedInit=transmittedInit)
        } 
        
        else {}
        
        oldstate=list(prop=prop)
        
        for(k in 1:length(prop)){
            oldstate=c(oldstate,list(sequ[k,]))
        }
        
    }
    
    return(transmission)
    ## output : list that contains for one host his ID, when he infects another host,
    ## how many virions are transmitted, in what proportions and what sequences
}




####################################################################
#####   SERIES OF SIMULATIONS
####################################################################

## function used for generating transmissions
transmit=function(lambda,virusmax,kinetics,checkinfection,pop){
    # checkinfection : bolean vector of host population (TRUE if healthy, FALSE if infected)
    # id : infectious host identifier
    susceptible=pop[checkinfection]
    duration=max(kinetics$time)-min(kinetics$time)
    contact.time=sample(1:nrow(kinetics),rpois(1,lambda*duration*length(susceptible)),
    replace=FALSE)
    contact.host=sample(susceptible,length(contact.time),replace=FALSE)
    proba.transmission=pmin(1,kinetics$virus[contact.time]/virusmax)
    success=(runif(length(contact.time))<proba.transmission)
    transmission.times=kinetics$time[contact.time[success]]
    transmission.hosts=contact.host[success]
    transmission=list(times=transmission.times,hosts=transmission.hosts)
    return(transmission)
    ## output : transmission times and respective target hosts
}

## function for generating an outbreak, simulating the multi-host compositions of sequences, and simulating the observed sequences
outbreak.chain=function(vpar,tpar,spar,fpar,opar,chainlength,directory,pop){
    ## vpar: parameters of viral kinetics
    ## tpar: parameters of transmission process
    ## spar: parameters for genetic sequence evolution
    ## fpar: fitness rule
    ## opar: parameters for observation
    ## directory: directory where output will be saved
    ## chainlength : length of the transmission chain
    ## pop : vector with identifier of every host in the considered population
  
    system(paste("rm -r",directory))
    system(paste("mkdir",directory))
    idmin=1000
    checkinfection=rep(TRUE,length(pop)) # initialization of checkinfection
    k=1 # k is the position of the host in the transmission chain
    id=1 # id of the host
    checkinfection[id]=FALSE # the 1st host is now infected
    STdyn=rbind(c(0,idmin+id,spar$init$time))
    ## list that describes the transmission chain, initialized to host n°1 (patient 0)
    vK=viralKinetics(vpar$init,vpar$s,vpar$d,vpar$beta,
        vpar$delta,vpar$p,vpar$c,vpar$timelag,vpar$duration)
    ## viral kinetics for host n°1
    transdyn=transmit(tpar$lambda,tpar$virusmax,kinetics=vK,checkinfection=checkinfection,pop=pop)
    ## determination of transmission times and respective target hosts for the initial host
    
    trans = withinHostData(init=spar$init,kinetics=vK,mu=spar$mu,gamma=spar$gamma,
    tolerance=spar$tolerance,fitnessRule=fpar,observationTime=sample(vK$time[vK$time>=opar$timemin
    & vK$time<=opar$timemax],1),samplesize=sample(opar$nsequmin:opar$nsequmax,1),directory=directory,
    identifier=idmin+id,id_chain=idmin+k,transdyn=transdyn,transmissionRate=tpar$rate)
    ## trans contains the transmission object of every transmission time
    ## transmission object : object returned by withinHostData
    
    while(k<chainlength & length(trans)>0 & k<=length(pop)){
    ## construction of the transmission chain while it is not completed, there still are susceptible individuals and trans isn't NULL  
        k=k+1
        ## infection = transmission to a new host = change of host
        nextinfectiontimes=as.numeric(lapply(trans,function(u) u$transmittedTime)) # times of transmission
        nextinfection=order(nextinfectiontimes)[1] # selection of the minimum transmission time
        transinit=trans[[nextinfection]] # selection of the transmission object for this time
        id=transinit$target # id of the target host for the selected transmission time
        if (checkinfection[id]==TRUE){
          checkinfection[id]=FALSE
          STdyn=rbind(STdyn,c(transinit$source,idmin+id,transinit$transmittedTime)) # addition of the new host to the transmission chain
          print(c(k,STdyn[nrow(STdyn),]))
          transold=trans
          trans=NULL
          
          for(i in (1:length(transold))[-nextinfection]){
              trans=c(trans,list(transold[[i]]))
          }
          ## withdrawal of the selected transmission object from trans
          
          vpar$init$V=transinit$transmittedVirus
          ## in the new host, initial V is transmittedVirus
          vK=viralKinetics(vpar$init,vpar$s,vpar$d,vpar$beta,
              vpar$delta,vpar$p,vpar$c,vpar$timelag,vpar$duration)
          ## viral kinetics in the new host
          vK$time=vK$time+transinit$transmittedTime
          ## time begins after transmission
          transdyn=transmit(tpar$lambda,tpar$virusmax,kinetics=vK,checkinfection=checkinfection,pop=pop)
          ## transmission times from current host to other individuals
          trans = c(trans,withinHostData(init=c(list(time=transinit$transmittedTime),transinit$transmittedInit),
          kinetics=vK,mu=spar$mu,gamma=spar$gamma,tolerance=spar$tolerance,fitnessRule=fpar,
          observationTime=sample(vK$time[vK$time>=vK$time[1]+opar$timemin & vK$time<=vK$time[1]+opar$timemax],1),
          samplesize=sample(opar$nsequmin:opar$nsequmax,1),directory=directory,id_chain=idmin+k,identifier=idmin+id,transdyn=transdyn,
          transmissionRate=tpar$rate))
        }
    }
    write.table(STdyn,paste(directory,"tree.txt",sep=""))
    return(STdyn)
    ## output : STdyn is the chain of transmission (describes who infected who and when)
}



###########################
#### WITHOUT FITNESS RULE
directory="data-simul-chain1/"
vpar1=list(init=list(TT=1e12,I=0,V=1e03),s=0,d=0,beta=5.0e-05,delta=2.0,p=1e-06,c=0.5,timelag=0.01,duration=50)
tpar1=list(lambda=0.002,virusmax=1e06,rate=1e-02,nsusceptible=1000)
initstate1=c(list(time=0,prop=c(0.7,rep(0.3/50,50)),sequmajor=rep(1,1e03)),sapply(1:50,function(u) list(sample(c(2,3,4,rep(1,1e03-3)),1e03,replace=TRUE))))
spar1=list(init=initstate1,mu=5e-06,gamma=c(0.001,1000),tolerance=1e-05)
fpar1=function(sequ){ return(1) }
opar1=list(timemin=0,timemax=3,nsequmin=30,nsequmax=70)
pop=1:1000
temp1=outbreak.chain(vpar=vpar1,tpar=tpar1,spar=spar1,fpar=fpar1,opar=opar1,chainlength=20,directory=directory,pop)
source("plotMST.r")

fpar2=function(sequ){
  M=length(sequ)
  return(1*(1+9*mean(sequ[1:round(0.1*M)]==2))/
           (1+9*mean(sequ[1:round(0.1*M)]==3)))
}

directory="data-simul-chain2/"
temp2=outbreak.chain(vpar=vpar1,tpar=tpar1,spar=spar1,fpar=fpar2,opar=opar1,chainlength=20,directory=directory,pop)
source("plotMST.r")


