
rm(list=ls())
nsimu=5
## load from data_ref file
# t=2
data_beta2=read.table("data_beta_t2.csv",header=TRUE,sep=",")
data_delta2=read.table("data_delta_t2.csv",header=TRUE,sep=",")
data_p2=read.table("data_p_t2.csv",header=TRUE,sep=",")
data_c2=read.table("data_c_t2.csv",header=TRUE,sep=",")
data_mu2=read.table("data_mu_t2.csv",header=TRUE,sep=",")
data_gamma12=read.table("data_gamma1_t2.csv",header=TRUE,sep=",")
data_gamma22=read.table("data_gamma2_t2.csv",header=TRUE,sep=",")
# t=2.5
data_beta25=read.table("data_beta_t25.csv",header=TRUE,sep=",")
data_delta25=read.table("data_delta_t25.csv",header=TRUE,sep=",")
data_p25=read.table("data_p_t25.csv",header=TRUE,sep=",")
data_c25=read.table("data_c_t25.csv",header=TRUE,sep=",")
data_mu25=read.table("data_mu_t25.csv",header=TRUE,sep=",")
data_gamma125=read.table("data_gamma1_t25.csv",header=TRUE,sep=",")
data_gamma225=read.table("data_gamma2_t25.csv",header=TRUE,sep=",")
# t=3
data_beta3=read.table("data_beta_t3.csv",header=TRUE,sep=",")
data_delta3=read.table("data_delta_t3.csv",header=TRUE,sep=",")
data_p3=read.table("data_p_t3.csv",header=TRUE,sep=",")
data_c3=read.table("data_c_t3.csv",header=TRUE,sep=",")
data_mu3=read.table("data_mu_t3.csv",header=TRUE,sep=",")
data_gamma13=read.table("data_gamma1_t3.csv",header=TRUE,sep=",")
data_gamma23=read.table("data_gamma2_t3.csv",header=TRUE,sep=",")
# t=15
data_beta15=read.table("data_beta_t15.csv",header=TRUE,sep=",")
data_delta15=read.table("data_delta_t15.csv",header=TRUE,sep=",")
data_p15=read.table("data_p_t15.csv",header=TRUE,sep=",")
data_c15=read.table("data_c_t15.csv",header=TRUE,sep=",")
data_mu15=read.table("data_mu_t15.csv",header=TRUE,sep=",")
data_gamma115=read.table("data_gamma1_t15.csv",header=TRUE,sep=",")
data_gamma215=read.table("data_gamma2_t15.csv",header=TRUE,sep=",")

Graph=function(data,name){
  par=data[,2]
  div=data[,3]
  plot(par,div,xlab="parameter",ylab="p-distance",main=name)
}


#########################
#### Time Comparison ####
#########################

# beta
par(mfrow=c(2,2))
Graph(data_beta2,"Beta, t=2")
Graph(data_beta25,"Beta, t=2.5")
Graph(data_beta3,"Beta, t=3")
Graph(data_beta15,"Beta, t=15")

# delta
par(mfrow=c(2,2))
Graph(data_delta2,"Delta, t=2")
Graph(data_delta25,"Delta, t=2.5")
Graph(data_delta3,"Delta, t=3")
Graph(data_delta15,"Delta, t=15")

# p
par(mfrow=c(2,2))
Graph(data_p2, "p, t=2")
Graph(data_p25, "p, t=2.5")
Graph(data_p3, "p, t=3")
Graph(data_p15, "p, t=15")

# c
par(mfrow=c(2,2))
Graph(data_c2,"c, t=2")
Graph(data_c25,"c, t=2.5")
Graph(data_c3,"c, t=3")
Graph(data_c15,"c, t=15")

# mu
par(mfrow=c(2,2))
Graph(data_mu2,"Mu, t=2")
Graph(data_mu25,"Mu, t=2.5")
Graph(data_mu3,"Mu, t=3")
Graph(data_mu15,"Mu, t=15")

# gamma1
par(mfrow=c(2,2))
Graph(data_gamma12,"Gamma1, t=2")
Graph(data_gamma125,"Gamma1, t=2.5")
Graph(data_gamma13,"Gamma1, t=3")
Graph(data_gamma115,"Gamma1, t=15")

# gamma2
par(mfrow=c(2,2))
Graph(data_gamma22,"Gamma2, t=2")
Graph(data_gamma225,"Gamma2, t=2.5")
Graph(data_gamma23,"Gamma2, t=3")
Graph(data_gamma215,"Gamma2, t=15")


###########################
#### Linear Regression ####
###########################

## predictions plot
# t=2
par(mfrow=c(3,3))
Reg_lm(data_beta2,"predict","Beta")
invdel2=1/data_delta2[,2]
newdata_d2=data.frame(data_delta2[,1],invdel2,data_delta2[,3])
Reg_lm(newdata_d2,"predict","Delta")
Reg_lm(data_p2,"predict","p")
Reg_lm(data_c2,"predict","c")
Reg_lm(data_mu2,"predict","mu")
lng12=log(data_gamma12[,2],base=exp(1))
newdata_g12=data.frame(data_gamma12[,1],lng12,data_gamma12[,3])
Reg_lm(newdata_g12,"predict","gamma1")
Reg_lm(data_gamma22,"predict","gamma2")

# t=2.5
par(mfrow=c(3,3))
Reg_lm(data_beta25,"predict","Beta")
invdel25=1/data_delta25[,2]
newdata_d25=data.frame(data_delta25[,1],invdel25,data_delta25[,3])
Reg_lm(newdata_d25,"predict","Delta")
Reg_lm(data_p25,"predict","p")
Reg_lm(data_c25,"predict","c")
Reg_lm(data_mu25,"predict","mu")
lng125=log(data_gamma125[,2],base=exp(1))
newdata_g125=data.frame(data_gamma125[,1],lng125,data_gamma125[,3])
Reg_lm(newdata_g125,"predict","gamma1")
Reg_lm(data_gamma225,"predict","gamma2")

# t=3
par(mfrow=c(3,3))
Reg_lm(data_beta3,"predict","Beta")
invdel3=1/data_delta3[,2]
newdata_d3=data.frame(data_delta3[,1],invdel3,data_delta3[,3])
Reg_lm(newdata_d3,"predict","Delta")
Reg_lm(data_p3,"predict","p")
Reg_lm(data_c3,"predict","c")
Reg_lm(data_mu3,"predict","mu")
lng13=log(data_gamma13[,2],base=exp(1))
newdata_g13=data.frame(data_gamma13[,1],lng13,data_gamma13[,3])
Reg_lm(newdata_g13,"predict","gamma1")
Reg_lm(data_gamma23,"predict","gamma2")

# t=15
par(mfrow=c(3,3))
Reg_lm(data_beta15,"predict","Beta")
invdel15=1/data_delta15[,2]
newdata_d15=data.frame(data_delta15[,1],invdel15,data_delta15[,3])
Reg_lm(data_delta15,"predict","Delta")
Reg_lm(data_p15,"predict","p")
Reg_lm(data_c15,"predict","c")
Reg_lm(data_mu15,"predict","mu")
lng115=log(data_gamma115[,2],base=exp(1))
newdata_g115=data.frame(data_gamma115[,1],lng12,data_gamma115[,3])
Reg_lm(newdata_g115,"predict","gamma1")
Reg_lm(data_gamma215,"predict","gamma2")


## hypothesis check
# t=2
par(mfrow=c(2,2))
Reg_lm(data_beta2,"check","Beta")
Reg_lm(data_delta2,"check","Delta")
Reg_lm(data_p2,"check","p")
Reg_lm(data_c2,"check","c")
Reg_lm(data_mu2,"check","mu")
Reg_lm(data_gamma12,"check","gamma1")
Reg_lm(data_gamma22,"check","gamma2")

# t=2.5
par(mfrow=c(2,2))
Reg_lm(data_beta25,"check","Beta")
Reg_lm(data_delta25,"check","Delta")
Reg_lm(data_p25,"check","p")
Reg_lm(data_c25,"check","c")
Reg_lm(data_mu25,"check","mu")
Reg_lm(data_gamma125,"check","gamma1")
Reg_lm(data_gamma225,"check","gamma2")

# t=3
par(mfrow=c(2,2))
Reg_lm(data_beta3,"check","Beta")
Reg_lm(data_delta3,"check","Delta")
Reg_lm(data_p3,"check","p")
Reg_lm(data_c3,"check","c")
Reg_lm(data_mu3,"check","mu")
Reg_lm(data_gamma13,"check","gamma1")
Reg_lm(data_gamma23,"check","gamma2")

# t=15
par(mfrow=c(2,2))
Reg_lm(data_beta15,"check","Beta")
Reg_lm(data_delta15,"check","Delta")
Reg_lm(data_p15,"check","p")
Reg_lm(data_c15,"check","c")
Reg_lm(data_mu15,"check","mu")
Reg_lm(data_gamma115,"check","gamma1")
Reg_lm(data_gamma215,"check","gamma2")


##########################
#### Spline Smoothing ####
##########################

# t=2
par(mfrow=c(3,3))
Reg_spl(data_beta2,"Beta")
Reg_spl(data_delta2,"Delta")
Reg_spl(data_p2,"p")
Reg_spl(data_c2,"c")
Reg_spl(data_mu2,"mu")
Reg_spl(data_gamma12,"gamma1")
Reg_spl(data_gamma22,"gamma2")

# t=2.5
par(mfrow=c(3,3))
Reg_spl(data_beta25,"Beta")
Reg_spl(data_delta25,"Delta")
Reg_spl(data_p25,"p")
Reg_spl(data_c25,"c")
Reg_spl(data_mu25,"mu")
Reg_spl(data_gamma125,"gamma1")
Reg_spl(data_gamma225,"gamma2")

# t=3
par(mfrow=c(3,3))
Reg_spl(data_beta3,"Beta")
Reg_spl(data_delta3,"Delta")
Reg_spl(data_p3,"p")
Reg_spl(data_c3,"c")
Reg_spl(data_mu3,"mu")
Reg_spl(data_gamma13,"gamma1")
Reg_spl(data_gamma23,"gamma2")

# t=15
par(mfrow=c(3,3))
Reg_spl(data_beta15,"Beta")
Reg_spl(data_delta15,"Delta")
Reg_spl(data_p15,"p")
Reg_spl(data_c15,"c")
Reg_spl(data_mu15,"mu")
Reg_spl(data_gamma115,"gamma1")
Reg_spl(data_gamma215,"gamma2")


###########################
#### Kernel Regression ####
###########################

# t=2
par(mfrow=c(3,3))
Reg_ker(data_beta2,"Beta")
Reg_ker(data_delta2,"Delta")
Reg_ker(data_p2,"p")
Reg_ker(data_c2,"c")
Reg_ker(data_mu2,"mu")
Reg_ker(data_gamma12,"gamma1")
Reg_ker(data_gamma22,"gamma2")

# t=2.5
par(mfrow=c(3,3))
Reg_ker(data_beta25,"Beta")
Reg_ker(data_delta25,"Delta")
Reg_ker(data_p25,"p")
Reg_ker(data_c25,"c")
Reg_ker(data_mu25,"mu")
Reg_ker(data_gamma125,"gamma1")
Reg_ker(data_gamma225,"gamma2")

# t=3
par(mfrow=c(3,3))
Reg_ker(data_beta3,"Beta")
Reg_ker(data_delta3,"Delta")
Reg_ker(data_p3,"p")
Reg_ker(data_c3,"c")
Reg_ker(data_mu3,"mu")
Reg_ker(data_gamma13,"gamma1")
Reg_ker(data_gamma23,"gamma2")

# t=15
par(mfrow=c(3,3))
Reg_ker(data_beta15,"Beta")
Reg_ker(data_delta15,"Delta")
Reg_ker(data_p15,"p")
Reg_ker(data_c15,"c")
Reg_ker(data_mu15,"mu")
Reg_ker(data_gamma115,"gamma1")
Reg_ker(data_gamma215,"gamma2")

##########################
#### Loess Regression ####
##########################

# Regression
# t=2
par(mfrow=c(3,3))
Reg_loess(data_beta2,"Beta","predict")
Reg_loess(data_delta2,"Delta","predict")
Reg_loess(data_p2,"p","predict")
Reg_loess(data_c2,"c","predict")
Reg_loess(data_mu2,"mu","predict")
Reg_loess(data_gamma12,"gamma1","predict")
Reg_loess(data_gamma22,"gamma2","predict")

# t=2.5
par(mfrow=c(3,3))
Reg_loess(data_beta25,"Beta","predict")
Reg_loess(data_delta25,"Delta","predict")
Reg_loess(data_p25,"p","predict")
Reg_loess(data_c25,"c","predict")
Reg_loess(data_mu25,"mu","predict")
Reg_loess(data_gamma125,"gamma1","predict")
Reg_loess(data_gamma225,"gamma2","predict")

# t=3
par(mfrow=c(3,3))
Reg_loess(data_beta3,"Beta","predict")
Reg_loess(data_delta3,"Delta","predict")
Reg_loess(data_p3,"p","predict")
Reg_loess(data_c3,"c","predict")
Reg_loess(data_mu3,"mu","predict")
Reg_loess(data_gamma13,"gamma1","predict")
Reg_loess(data_gamma23,"gamma2","predict")

# t=15
par(mfrow=c(3,3))
Reg_loess(data_beta15,"Beta","predict")
Reg_loess(data_delta15,"Delta","predict")
Reg_loess(data_p15,"p","predict")
Reg_loess(data_c15,"c","predict")
Reg_loess(data_mu15,"mu","predict")
Reg_loess(data_gamma115,"gamma1","predict")
Reg_loess(data_gamma215,"gamma2","predict")

# Hypothesis check
# t=2
par(mfrow=c(2,2))
l1=Reg_loess(data_beta2,"Beta","check")
par(mfrow=c(2,2))
l2=Reg_loess(data_delta2,"Delta","check")
par(mfrow=c(2,2))
l3=Reg_loess(data_p2,"p","check")
par(mfrow=c(2,2))
l4=Reg_loess(data_c2,"c","check")
par(mfrow=c(2,2))
l5=Reg_loess(data_mu2,"mu","check")
par(mfrow=c(2,2))
l6=Reg_loess(data_gamma12,"gamma1","check")
par(mfrow=c(2,2))
l7=Reg_loess(data_gamma22,"gamma2","check")
S1=sum(l1$p.value,l2$p.value,l3$p.value,l4$p.value,l5$p.value,l6$p.value,l7$p.value)

# t=2.5
par(mfrow=c(2,2))
m1=Reg_loess(data_beta25,"Beta","check")
par(mfrow=c(2,2))
m2=Reg_loess(data_delta25,"Delta","check")
par(mfrow=c(2,2))
m3=Reg_loess(data_p25,"p","check")
par(mfrow=c(2,2))
m4=Reg_loess(data_c25,"c","check")
par(mfrow=c(2,2))
m5=Reg_loess(data_mu25,"mu","check")
par(mfrow=c(2,2))
m6=Reg_loess(data_gamma125,"gamma1","check")
par(mfrow=c(2,2))
m7=Reg_loess(data_gamma225,"gamma2","check")
S2=sum(m1$p.value,m2$p.value,m3$p.value,m4$p.value,m5$p.value,m6$p.value,m7$p.value)

# t=3
par(mfrow=c(2,2))
n1=Reg_loess(data_beta3,"Beta","check")
par(mfrow=c(2,2))
n2=Reg_loess(data_delta3,"Delta","check")
par(mfrow=c(2,2))
n3=Reg_loess(data_p3,"p","check")
par(mfrow=c(2,2))
n4=Reg_loess(data_c3,"c","check")
par(mfrow=c(2,2))
n5=Reg_loess(data_mu3,"mu","check")
par(mfrow=c(2,2))
n6=Reg_loess(data_gamma13,"gamma1","check")
par(mfrow=c(2,2))
n7=Reg_loess(data_gamma23,"gamma2","check")
S3=sum(n1$p.value,n2$p.value,n3$p.value,n4$p.value,n5$p.value,n6$p.value,n7$p.value)

# t=15
par(mfrow=c(2,2))
p1=Reg_loess(data_beta15,"Beta","check")
par(mfrow=c(2,2))
p2=Reg_loess(data_delta15,"Delta","check")
par(mfrow=c(2,2))
p3=Reg_loess(data_p15,"p","check")
par(mfrow=c(2,2))
p4=Reg_loess(data_c15,"c","check")
par(mfrow=c(2,2))
p5=Reg_loess(data_mu15,"mu","check")
par(mfrow=c(2,2))
p6=Reg_loess(data_gamma115,"gamma1","check")
par(mfrow=c(2,2))
p7=Reg_loess(data_gamma215,"gamma2","check")
S4=sum(p1$p.value,p2$p.value,p3$p.value,p4$p.value,p5$p.value,p6$p.value,p7$p.value)

########################
#### Gam Regression ####
########################

# t=2
par(mfrow=c(3,3))
Reg_gam(data_beta2,"Beta")
Reg_gam(data_delta2,"Delta")
Reg_gam(data_p2,"p")
Reg_gam(data_c2,"c")
Reg_gam(data_mu2,"mu")
Reg_gam(data_gamma12,"gamma1")
Reg_gam(data_gamma22,"gamma2")

# t=2.5
par(mfrow=c(3,3))
Reg_gam(data_beta25,"Beta")
Reg_gam(data_delta25,"Delta")
Reg_gam(data_p25,"p")
Reg_gam(data_c25,"c")
Reg_gam(data_mu25,"mu")
Reg_gam(data_gamma125,"gamma1")
Reg_gam(data_gamma225,"gamma2")

# t=3
par(mfrow=c(3,3))
Reg_gam(data_beta3,"Beta")
Reg_gam(data_delta3,"Delta")
Reg_gam(data_p3,"p")
Reg_gam(data_c3,"c")
Reg_gam(data_mu3,"mu")
Reg_gam(data_gamma13,"gamma1")
Reg_gam(data_gamma23,"gamma2")

# t=15
par(mfrow=c(3,3))
Reg_gam(data_beta15,"Beta")
Reg_gam(data_delta15,"Delta")
Reg_gam(data_p15,"p")
Reg_gam(data_c15,"c")
Reg_gam(data_mu15,"mu")
Reg_gam(data_gamma115,"gamma1")
Reg_gam(data_gamma215,"gamma2")




