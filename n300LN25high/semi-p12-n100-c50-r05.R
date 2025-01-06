#Simulation Study for Right-censored Left-truncated Data using shared frailty and Cox multiplicative hazards 
#functions with the baseline hazards functions estimated using B-spline method for 3 transitions:
#Transition 1: stage 1 ---> Non-terminal event
#Transition 2: stage 1 ---> Terminal event
#Transition 3: Non-terminal event---> Terminal event 

#######Run the preamble  ###################################################
#For in ARC running:
# getwd()
# setwd("")
rm(list=ls())
#For on computer running:
setwd(getwd())
#replace this directory with the one of the unzipped files
# setwd("/Users/fatemehmahmoudi/Desktop/ARC-BSpline/ARC-BS-BAR")
#setwd("/Users/fmahmoudi/Documents/Research/GNST2023/Copy-Codes-300-rep/n300LN25low")
source("sim-functions-bp-revised.R")
# source("varsel-functions-bp-revised-2.R")
source("varsel-bp-2.R")

###########################################################################
#For running bactehes on ARC:
#####This part is for running batches in ARC:#######
args = commandArgs(trailingOnly = TRUE)
start_seed = as.numeric(args[1])
end_seed = as.numeric(args[2])
seeds = seq(start_seed, end_seed)
cat("Running simulations for seeds:", seeds, "\n")
###################n=?,rho=?################
# options(digits=3,width=200)
# 10 variables considered for each transition with the following values:
#seeds=c(1:50)
# seeds=c(21:40)
# seeds=c(41:60)
# seeds=c(61:80)
# seeds=c(81:100)


#Sample size:
n=300.   #We test n=100 with pn=12, 300 with pn=15 and 500 with pn=16
#censoring="low"
censoring="high"

#If we want to make GM comparable to Gamma distribution setting used in the manuscript, we need to find sigma^2's corresponding to 
#theta values equal to 0.25v and 1. 
#If we calculate those, corresponding to theta=0.25, we will have sigma^2=0.188.
#corresponding to theta=1, we will have sigma^2=0.481.

#theta.given=0.481 #comparable to theta=1.00.
#we generate LN with E(X)=1 and variance =0.2231 corresponding to theta=0.25 and 0.6931 corresponding to theta=1:
theta.given=sqrt(0.2231) #comparable to theta=0.25. theta.given is sigma^2 in LN.
# theta.given=sqrt(0.6931)
# frailty_distribution="Gamma"
frailty_distribution="LogNormal"
lt.type="unif"




if (n==100){
  beta1.true = c(-0.8,1,1,0.9,rep(0,8))
  beta2.true = c(1,1,1,0.9,rep(0,8))
  beta3.true = c(-1,1,0.9,1,rep(0,8))
}
if (n==300){
  beta1.true = c(-0.8,1,1,0.9,rep(0,11))
  beta2.true = c(1,1,1,0.9,rep(0,11))
  beta3.true = c(-1,1,0.9,1,rep(0,11))
}
if (n==500){
  beta1.true = c(-0.8,1,1,0.9,rep(0,12))
  beta2.true = c(1,1,1,0.9,rep(0,12))
  beta3.true = c(-1,1,0.9,1,rep(0,12))
}

weibull.param.log=c(-4,0.18,-4,0.2,-11,1.7,-1.4)


if (censoring=="low"){
  c1=1000
  c2=1000
}
if (censoring=="high"){
  c1=45
  c2=45
}

frailty=TRUE
# risks.status="overlapping"
risks.status="fullyshared"
#grouping effect study or individual variable selection:
grp.or.indiv="indiv"

method="optim" #this one is preferred. Some settings work with "nlm" only, though.
# method="nlm"
yvalslength=500 # #of some y values that are required to be generated while approximating baseline hazards using bernstein polynomials. 
lt1=0
lt2=0.3

m=c(2,2,3) #degree of Bernstein polynomial for three transitions corresponding to 3 transitions 1,2,and 3.
myseed=paste(min(seeds),max(seeds),sep=":")

if (grp.or.indiv=="indiv"){
  # rho1=0.2
  rho2=0.5
  # rho3=0.8
}

if (grp.or.indiv=="grp"){
  rho1=0.8
  rho2=0.9
  rho3=0.95
}

#group effect study parameters:
G1=1:2
G2=3:4
G3=5:7
G4=8:10
G1.distr="normal";G2.distr="Bernoulli";G3.distr="normal";G4.distr="Bernoulli"
mu.bernoulli=0.5

#not required:
# lambdavec.lasso=seq(8.9,9.4,0.05)
# lambdavec.ada=seq(8.9,9.4,0.05)
# lambdavec.bar=c(1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6)



if (n==100){
  lambdavec.lasso=seq(8,9,0.15)
  lambdavec.ada=seq(8,9,0.15)
  lambdavec.bar=seq(0.5,1.5,0.15)
}

if (n==300){
  lambdavec.lasso=seq(8.763,9.763,0.15)
  lambdavec.ada=seq(8.763,9.763,0.15)
  lambdavec.bar=seq(0.612,1.612,0.15)
}

if (n==500){
  lambdavec.lasso=seq(8.984,9.984,0.15)
  lambdavec.ada=seq(8.984,9.984,0.15)
  lambdavec.bar=seq(0.644,1.644,0.15)
}



#indiv:
# lambdavec.lasso=seq(5,6.5,0.25)
# lambdavec.ada=seq(2,3.5,0.25)
# lambdavec.bar=seq(1.5,3,0.25)
# n100p12:
# lambdavec.lasso=seq(8,9,0.15)
# lambdavec.ada=seq(8,9,0.15)
# lambdavec.bar=seq(0.5,1.5,0.15)
# n300p15:
# lambdavec.lasso=seq(8.763,9.763,0.15)
# lambdavec.ada=seq(8.763,9.763,0.15)
# # lambdavec.lasso=seq(5,6.5,1)
# # lambdavec.ada=seq(2,3.5,1)
# # lambdavec.bar=seq(1,2,0.25)
# # lambdavec.bar=seq(1,2.5,0.25)
# lambdavec.bar=seq(0.612,1.612,0.15)
#grp:
# lambdavec.lasso=seq(5,6,0.25)
# lambdavec.ada=seq(4,5,0.25)
# lambdavec.bar=seq(0.75,1.75,0.25)

xi=10
tol=1e-04#criterion for selecting zero and nonzero in BAR.
optim_starting_values=matrix(0,length(seeds),11)

E_ZZprime.for.oracle=list()
censoring.observ.rates=matrix(0,length(seeds),4)#first col: both censored, 2nd col:nonterminal not observed, 3rd col: terminal censored, 4th col: both observed.
MSE.calculator.oracle=vector()
p=length(beta1.true)
dimension=vector()

G1.percentage.bar=vector()
G2.percentage.bar=vector()
G3.percentage.bar=vector()
G4.percentage.bar=vector()

G1.percentage.lasso=vector()
G2.percentage.lasso=vector()
G3.percentage.lasso=vector()
G4.percentage.lasso=vector()

G1.percentage.ada=vector()
G2.percentage.ada=vector()
G3.percentage.ada=vector()
G4.percentage.ada=vector()

G1.for.that.lam=vector()
G2.for.that.lam=vector()
G3.for.that.lam=vector()
G4.for.that.lam=vector()

BAR.updated.theta=matrix(0,11,length(seeds))
TP.bar=vector()
FP.bar=vector()
FN.bar=vector()
MSE.bar=vector()
GCV.bar=vector()
theta.bar=matrix(0,11,length(seeds))
GCV.selected.bar=matrix(0,length(lambdavec.bar),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.bar))
FP.for.that.lam=vector()
MSE.for.that.lam=vector()
theta.for.that.lam.bar=matrix(0,11,length(lambdavec.bar))
betavarsel.bar=matrix(0,3*p,length(lambdavec.bar))
beta.selected.bar=matrix(0,3*p,length(seeds))

adalasso.updated.theta=matrix(0,11,length(seeds))
TP.ada=vector()
FP.ada=vector()
FN.ada=vector()
MSE.ada=vector()
theta.ada=matrix(0,11,length(seeds))
GCV.AdaLasso=vector()
GCV.selected.ada=matrix(0,length(lambdavec.ada),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.ada))
FP.for.that.lam=vector()
FN.for.that.lam=vector()
MSE.for.that.lam=vector()
theta.for.that.lam.AdaLasso=matrix(0,11,length(lambdavec.ada))
betavarsel.AdaLasso=matrix(0,3*p,length(lambdavec.ada))
beta.selected.ada=matrix(0,3*p,length(seeds))


lasso.updated.theta=matrix(0,11,length(seeds))
TP.lasso=vector()
FP.lasso=vector()
FN.lasso=vector()
MSE.lasso=vector()
GCV.lasso=vector()
theta.lasso=matrix(0,11,length(seeds))
GCV.selected.lasso=matrix(0,length(lambdavec.lasso),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.lasso))
FP.for.that.lam=vector()
FN.for.that.lam=vector()
MSE.for.that.lam=vector()
theta.for.that.lam.lasso=matrix(0,11,length(lambdavec.lasso))
betavarsel.lasso=matrix(0,3*p,length(lambdavec.lasso))
beta.selected.lasso=matrix(0,3*p,length(seeds))

#Update theta (nuisance params with each peanlty):
starttime=date()
d=repeat.get.last.result.bp.update.theta(seeds,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho2,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,lambdavec.bar,lambdavec.lasso,lambdavec.ada,theta.given,frailty_distribution)
endtime=date()

#Not uopdating nuisance params with penalties (the basic version before May 2024 update after revision from JNS):
# starttime=date()
# d=repeat.get.last.result.bp(seeds,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho2,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,lambdavec.bar,lambdavec.lasso,lambdavec.ada,theta.given,frailty_distribution)
# endtime=date()

filename=filename="semi-p12-n100-c50-sg-r05"

sink(filename)
cat("Simulation results for the Param-Bernstein Poly. Semi-parametric Semi-competing model:", "\n")

cat("sample size              =",n,"\n",
    "censoring value          =",c1,"\n",
    "number of simulations    =",length(seeds),"\n",
    "pn                       =",length(beta1.true),"\n",
    "grporindiv               =",grp.or.indiv,"\n",
    "strongorweak             =",sum(beta1.true),"\n\n\n\n\n")

cat("Simulation results for nuisance parameters based on unpenalized likelihood in BP model:", "\n\n")

#UNPENALIZED ESTIMATE NUISANCE PARAMS:
#on LOG scale:
cat("Estimation results based on the log scale for unpenalized way:", "\n \n")

cat(" True values of log of tau_1                  =",weibull.param.log[1],"\n",
    "True values of log of alpha_1                 =",weibull.param.log[2],"\n",
    "True values of log of tau_2                   =",weibull.param.log[3],"\n",
    "True values of log of alpha_2                 =",weibull.param.log[4],"\n",
    "True values of log of tau_3                   =",weibull.param.log[5],"\n",
    "True values of log of alpha_3                 =",weibull.param.log[6],"\n",
    "True values of log of theta                   =",log(theta.given),"\n\n\n\n"
)

# cat("Average starting value in optim function in log-scale:", "\n \n")
# cat("starting value in optim function for tau_1 in log scale                               =",mean(d$optim_starting_values[,1]),"\n",
#     "starting value in optim function for alpha_1 in log scale                             =",mean(d$optim_starting_values[,2]),"\n",
#     "starting value in optim function for tau_2 in log scale                               =",mean(d$optim_starting_values[,3]),"\n",
#     "starting value in optim function for alpha_2 in log scale                             =",mean(d$optim_starting_values[,4]),"\n",
#     "starting value in optim function for tau_3 in log scale                               =",mean(d$optim_starting_values[,5]),"\n",
#     "starting value in optim function for alpha_3 in log scale                             =",mean(d$optim_starting_values[,6]),"\n",
#     "starting value in optim function for theta (frailty para) in log scale                =",mean(d$logtheta.init),"\n\n\n")
# 

cat("Mean Estimate of log of phi_10 across all the simulations                               =",d$theta.bar.mean[1],"\n",
    "Mean Estimate of log of phi_11 across all the simulations                               =",d$theta.bar.mean[2],"\n",
    "Mean Estimate of log of phi_12 across all the simulations                               =",d$theta.bar.mean[3],"\n",
    "Mean Estimate of log of phi_20 across all the simulations                               =",d$theta.bar.mean[4],"\n",
    "Mean Estimate of log of phi_21 across all the simulations                               =",d$theta.bar.mean[5],"\n",
    "Mean Estimate of log of phi_22 across all the simulations                               =",d$theta.bar.mean[6],"\n",
    "Mean Estimate of log of phi_30 across all the simulations                               =",d$theta.bar.mean[7],"\n",
    "Mean Estimate of log of phi_31 across all the simulations                               =",d$theta.bar.mean[8],"\n",
    "Mean Estimate of log of phi_32 across all the simulations                               =",d$theta.bar.mean[9],"\n",
    "Mean Estimate of log of phi_33 across all the simulations                               =",d$theta.bar.mean[10],"\n",
    "Mean Estimate of log of frailty pramaeter across all the simulations                    =",d$theta.bar.mean[11],"\n\n\n\n"
    
)

cat("Standard deviation of log of phi_10 across all the simulations                               =",sd(d$theta.bar[1,]),"\n",
    "Standard deviation of log of phi_11 across all the simulations                               =",sd(d$theta.bar[2,]),"\n",
    "Standard deviation of log of phi_12 across all the simulations                               =",sd(d$theta.bar[3,]),"\n",
    "Standard deviation of log of phi_20 across all the simulations                               =",sd(d$theta.bar[4,]),"\n",
    "Standard deviation of log of phi_21 across all the simulations                               =",sd(d$theta.bar[5,]),"\n",
    "Standard deviation of log of phi_22 across all the simulations                               =",sd(d$theta.bar[6,]),"\n",
    "Standard deviation of log of phi_30  across all the simulations                              =",sd(d$theta.bar[7,]),"\n",
    "Standard deviation of log of phi_31  across all the simulations                              =",sd(d$theta.bar[8,]),"\n",
    "Standard deviation of log of phi_32  across all the simulations                              =",sd(d$theta.bar[9,]),"\n",
    "Standard deviation of log of phi_33  across all the simulations                              =",sd(d$theta.bar[10,]),"\n",
    "Standard deviation of log of frailty pramaetr  across all the simulations                    =",sd(d$theta.bar[11,]),"\n\n\n"
    
)


#on ORIGINAL scale:
cat("Estimation results based on the original scale:", "\n \n")
cat(" True values of log of tau_1                   =",exp(weibull.param.log[1]),"\n",
    "True values of log of alpha_1                 =",exp(weibull.param.log[2]),"\n",
    "True values of log of tau_2                   =",exp(weibull.param.log[3]),"\n",
    "True values of log of alpha_2                 =",exp(weibull.param.log[4]),"\n",
    "True values of log of tau_3                   =",exp(weibull.param.log[5]),"\n",
    "True values of log of alpha_3                 =",exp(weibull.param.log[6]),"\n",
    "True values of log of theta                   =",theta.given,"\n\n\n\n"
)
cat("Average starting value in optim function in original-scale:", "\n \n")
cat("starting value in optim function for phi_10 in original scale                               =",mean(exp(d$optim_starting_values[,1])),"\n",
    "starting value in optim function for phi_11 in original scale                               =",mean(exp(d$optim_starting_values[,2])),"\n",
    "starting value in optim function for phi_12 in original scale                               =",mean(exp(d$optim_starting_values[,3])),"\n",
    "starting value in optim function for phi_20 in original scale                               =",mean(exp(d$optim_starting_values[,4])),"\n",
    "starting value in optim function for phi_21 in original scale                               =",mean(exp(d$optim_starting_values[,5])),"\n",
    "starting value in optim function for phi_22 in original scale                               =",mean(exp(d$optim_starting_values[,6])),"\n",
    "starting value in optim function for phi_30 in original scale                               =",mean(exp(d$optim_starting_values[,7])),"\n",
    "starting value in optim function for phi_31 in original scale                               =",mean(exp(d$optim_starting_values[,8])),"\n",
    "starting value in optim function for phi_32 in original scale                               =",mean(exp(d$optim_starting_values[,9])),"\n",
    "starting value in optim function for phi_33 in original scale                               =",mean(exp(d$optim_starting_values[,10])),"\n",
    "starting value in optim function for frailty param in original scale                        =",mean(exp(d$optim_starting_values[,11])),"\n\n\n\n"
    )


cat("Mean Estimate of of phi_10 across all the simulations                               =",mean(exp(d$theta.bar[1,])),"\n",
    "Mean Estimate of of phi_11 across all the simulations                               =",mean(exp(d$theta.bar[2,])),"\n",
    "Mean Estimate of of phi_12 across all the simulations                               =",mean(exp(d$theta.bar[3,])),"\n",
    "Mean Estimate of of phi_20 across all the simulations                               =",mean(exp(d$theta.bar[4,])),"\n",
    "Mean Estimate of of phi_21 across all the simulations                               =",mean(exp(d$theta.bar[5,])),"\n",
    "Mean Estimate of of phi_22 across all the simulations                               =",mean(exp(d$theta.bar[6,])),"\n",
    "Mean Estimate of of phi_30 across all the simulations                               =",mean(exp(d$theta.bar[7,])),"\n",
    "Mean Estimate of of phi_31 across all the simulations                               =",mean(exp(d$theta.bar[8,])),"\n",
    "Mean Estimate of of phi_32 across all the simulations                               =",mean(exp(d$theta.bar[9,])),"\n",
    "Mean Estimate of of phi_33 across all the simulations                               =",mean(exp(d$theta.bar[10,])),"\n",
    "Mean Estimate of of frailty param across all the simulations                        =",mean(exp(d$theta.bar[11,])),"\n\n\n\n"
    
    
)

cat("Standard deviation of of phi_10 across all the simulations                                =",sd(exp(d$theta.bar[1,])),"\n",
    "Standard deviation of of phi_11 across all the simulations                                =",sd(exp(d$theta.bar[2,])),"\n",
    "Standard deviation of of phi_12 across all the simulations                                =",sd(exp(d$theta.bar[3,])),"\n",
    "Standard deviation of of phi_20 across all the simulations                                =",sd(exp(d$theta.bar[4,])),"\n",
    "Standard deviation of of phi_21 across all the simulations                                =",sd(exp(d$theta.bar[5,])),"\n",
    "Standard deviation of of phi_22 across all the simulations                                =",sd(exp(d$theta.bar[6,])),"\n",
    "Standard deviation of of phi_30 across all the simulations                                =",sd(exp(d$theta.bar[7,])),"\n", 
    "Standard deviation of of phi_31 across all the simulations                                =",sd(exp(d$theta.bar[8,])),"\n",
    "Standard deviation of of phi_32 across all the simulations                                =",sd(exp(d$theta.bar[9,])),"\n",
    "Standard deviation of of phi_33 across all the simulations                                =",sd(exp(d$theta.bar[10,])),"\n",
    "Standard deviation of of frailty param across all the simulations                         =",sd(exp(d$theta.bar[11,])),"\n\n\n\n"
    )
cat("Bias Estimate of theta (frailty parameter) across all simulations                 =",mean(exp(d$theta.bar[11,]))-theta.given,"\n\n\n"
)



#theta vector containing three weibull parameters and one frailty parameter is updated after variable selection
#using eacch penalty.

#UPDATED ESTIMATES of  NUISANCE PARAMS BASED on ****BAR:
#on LOG scale:
cat("Estimation results using BAR Penalty based on the log scale:", "\n \n")


cat("BAR Updated Mean Estimate of log of phi_10 across all the simulations                               =",mean(d$BAR.updated.theta[1,]),"\n",
    "BAR Updated Mean Estimate of log of phi_11 across all the simulations                               =",mean(d$BAR.updated.theta[2,]),"\n",
    "BAR Updated Mean Estimate of log of phi_12 across all the simulations                               =",mean(d$BAR.updated.theta[3,]),"\n",
    "BAR Updated Mean Estimate of log of phi_20 across all the simulations                               =",mean(d$BAR.updated.theta[4,]),"\n",
    "BAR Updated Mean Estimate of log of phi_21 across all the simulations                               =",mean(d$BAR.updated.theta[5,]),"\n",
    "BAR Updated Mean Estimate of log of phi_22 across all the simulations                               =",mean(d$BAR.updated.theta[6,]),"\n",
    "BAR Updated Mean Estimate of log of phi_30 across all the simulations                               =",mean(d$BAR.updated.theta[7,]),"\n",
    "BAR Updated Mean Estimate of log of phi_31 across all the simulations                               =",mean(d$BAR.updated.theta[8,]),"\n",
    "BAR Updated Mean Estimate of log of phi_32 across all the simulations                               =",mean(d$BAR.updated.theta[9,]),"\n",
    "BAR Updated Mean Estimate of log of phi_33 across all the simulations                               =",mean(d$BAR.updated.theta[10,]),"\n",
    "BAR Updated Mean Estimate of log of frailty param across all the simulations                        =",mean(d$BAR.updated.theta[11,]),"\n\n\n"
    
)


cat("BAR Updated Standard deviation of log of phi_10 across all the simulations                               =",sd(d$BAR.updated.theta[1,]),"\n",
    "BAR Updated Standard deviation of log of phi_11 across all the simulations                               =",sd(d$BAR.updated.theta[2,]),"\n",
    "BAR Updated Standard deviation of log of phi_12 across all the simulations                               =",sd(d$BAR.updated.theta[3,]),"\n",
    "BAR Updated Standard deviation of log of phi_20 across all the simulations                               =",sd(d$BAR.updated.theta[4,]),"\n",
    "BAR Updated Standard deviation of log of phi_21 across all the simulations                               =",sd(d$BAR.updated.theta[5,]),"\n",
    "BAR Updated Standard deviation of log of phi_22 across all the simulations                               =",sd(d$BAR.updated.theta[6,]),"\n",
    "BAR Updated Standard deviation of log of phi_30 across all the simulations                               =",sd(d$BAR.updated.theta[7,]),"\n",
    "BAR Updated Standard deviation of log of phi_31 across all the simulations                               =",sd(d$BAR.updated.theta[8,]),"\n",
    "BAR Updated Standard deviation of log of phi_32 across all the simulations                               =",sd(d$BAR.updated.theta[9,]),"\n",
    "BAR Updated Standard deviation of log of phi_33 across all the simulations                               =",sd(d$BAR.updated.theta[10,]),"\n",
    "BAR Updated Standard deviation of log of frailty param across all the simulations                        =",sd(d$BAR.updated.theta[11,]),"\n\n\n\n"

)

cat("BAR Updated Bias Estimate of log of theta (frailty parameter) across all the simulations             =",mean(d$BAR.updated.theta[11,])-log(theta.given),"\n\n\n"
)

#on ORIGINAL scale:


cat("BAR Updated Mean Estimate of of phi_10 across all the simulations                               =",mean(exp(d$BAR.updated.theta[1,])),"\n",
    "BAR Updated Mean Estimate of of phi_11 across all the simulations                               =",mean(exp(d$BAR.updated.theta[2,])),"\n",
    "BAR Updated Mean Estimate of of phi_12 across all the simulations                               =",mean(exp(d$BAR.updated.theta[3,])),"\n",
    "BAR Updated Mean Estimate of of phi_20 across all the simulations                               =",mean(exp(d$BAR.updated.theta[4,])),"\n",
    "BAR Updated Mean Estimate of of phi_21 across all the simulations                               =",mean(exp(d$BAR.updated.theta[5,])),"\n",
    "BAR Updated Mean Estimate of of phi_22 across all the simulations                               =",mean(exp(d$BAR.updated.theta[6,])),"\n",
    "BAR Updated Mean Estimate of of phi_30 across all the simulations                               =",mean(exp(d$BAR.updated.theta[7,])),"\n",
    "BAR Updated Mean Estimate of of phi_31 across all the simulations                               =",mean(exp(d$BAR.updated.theta[8,])),"\n",
    "BAR Updated Mean Estimate of of phi_32 across all the simulations                               =",mean(exp(d$BAR.updated.theta[9,])),"\n",
    "BAR Updated Mean Estimate of of phi_33 across all the simulations                               =",mean(exp(d$BAR.updated.theta[10,])),"\n",
    "BAR Updated Mean Estimate of of frailty param across all the simulations                        =",mean(exp(d$BAR.updated.theta[11,])),"\n\n\n\n"
    
)

cat("BAR Updated Standard deviation of of phi_10 across all the simulations                               =",sd(exp(d$BAR.updated.theta[1,])),"\n",
    "BAR Updated Standard deviation of of phi_11 across all the simulations                               =",sd(exp(d$BAR.updated.theta[2,])),"\n",
    "BAR Updated Standard deviation of of phi_12 across all the simulations                               =",sd(exp(d$BAR.updated.theta[3,])),"\n",
    "BAR Updated Standard deviation of of phi_20 across all the simulations                               =",sd(exp(d$BAR.updated.theta[4,])),"\n",
    "BAR Updated Standard deviation of of phi_21 across all the simulations                               =",sd(exp(d$BAR.updated.theta[5,])),"\n",
    "BAR Updated Standard deviation of of phi_22 across all the simulations                               =",sd(exp(d$BAR.updated.theta[6,])),"\n",
    "BAR Updated Standard deviation of of phi_30 across all the simulations                               =",sd(exp(d$BAR.updated.theta[7,])),"\n",
    "BAR Updated Standard deviation of of phi_31 across all the simulations                               =",sd(exp(d$BAR.updated.theta[8,])),"\n",
    "BAR Updated Standard deviation of of phi_32 across all the simulations                               =",sd(exp(d$BAR.updated.theta[9,])),"\n",
    "BAR Updated Standard deviation of of phi_33 across all the simulations                               =",sd(exp(d$BAR.updated.theta[10,])),"\n",
    "BAR Updated Standard deviation of of frailty param across all the simulations                        =",sd(exp(d$BAR.updated.theta[11,])),"\n\n\n\n"
    
)


cat("BAR Updated Bias Estimate of theta (frailty parameter) across all simulations                 =",mean(exp(d$BAR.updated.theta[11,]))-theta.given,"\n\n\n"
)





#UPDATED ESTIMATES of  NUISANCE PARAMS BASED on ****LASSO:
#on LOG scale:
cat("Estimation results using LASSO Penalty based on the log scale:", "\n \n")





cat("LASSO Updated Mean Estimate of log of phi_10 across all the simulations                               =",mean(d$lasso.updated.theta[1,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_11 across all the simulations                               =",mean(d$lasso.updated.theta[2,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_12 across all the simulations                               =",mean(d$lasso.updated.theta[3,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_20 across all the simulations                               =",mean(d$lasso.updated.theta[4,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_21 across all the simulations                               =",mean(d$lasso.updated.theta[5,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_22 across all the simulations                               =",mean(d$lasso.updated.theta[6,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_30 across all the simulations                               =",mean(d$lasso.updated.theta[7,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_31 across all the simulations                               =",mean(d$lasso.updated.theta[8,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_32 across all the simulations                               =",mean(d$lasso.updated.theta[9,]),"\n",
    "LASSO Updated Mean Estimate of log of phi_33 across all the simulations                               =",mean(d$lasso.updated.theta[10,]),"\n",
    "LASSO Updated Mean Estimate of log of frailty param across all the simulations                        =",mean(d$lasso.updated.theta[11,]),"\n\n\n"
    
)


cat("LASSO Updated Standard deviation of log of phi_10 across all the simulations                               =",sd(d$lasso.updated.theta[1,]),"\n",
    "LASSO Updated Standard deviation of log of phi_11 across all the simulations                               =",sd(d$lasso.updated.theta[2,]),"\n",
    "LASSO Updated Standard deviation of log of phi_12 across all the simulations                               =",sd(d$lasso.updated.theta[3,]),"\n",
    "LASSO Updated Standard deviation of log of phi_20 across all the simulations                               =",sd(d$lasso.updated.theta[4,]),"\n",
    "LASSO Updated Standard deviation of log of phi_21 across all the simulations                               =",sd(d$lasso.updated.theta[5,]),"\n",
    "LASSO Updated Standard deviation of log of phi_22 across all the simulations                               =",sd(d$lasso.updated.theta[6,]),"\n",
    "LASSO Updated Standard deviation of log of phi_30 across all the simulations                               =",sd(d$lasso.updated.theta[7,]),"\n",
    "LASSO Updated Standard deviation of log of phi_31 across all the simulations                               =",sd(d$lasso.updated.theta[8,]),"\n",
    "LASSO Updated Standard deviation of log of phi_32 across all the simulations                               =",sd(d$lasso.updated.theta[9,]),"\n",
    "LASSO Updated Standard deviation of log of phi_33 across all the simulations                               =",sd(d$lasso.updated.theta[10,]),"\n",
    "LASSO Updated Standard deviation of log of frailty param across all the simulations                        =",sd(d$lasso.updated.theta[11,]),"\n\n\n\n"
    
)

cat("LASSO Updated Bias Estimate of log of theta (frailty parameter) across all the simulations             =",mean(d$lasso.updated.theta[11,])-log(theta.given),"\n\n\n"
)


#on ORIGINAL scale:


cat("LASSO Updated Mean Estimate of of phi_10 across all the simulations                               =",mean(exp(d$lasso.updated.theta[1,])),"\n",
    "LASSO Updated Mean Estimate of of phi_11 across all the simulations                               =",mean(exp(d$lasso.updated.theta[2,])),"\n",
    "LASSO Updated Mean Estimate of of phi_12 across all the simulations                               =",mean(exp(d$lasso.updated.theta[3,])),"\n",
    "LASSO Updated Mean Estimate of of phi_20 across all the simulations                               =",mean(exp(d$lasso.updated.theta[4,])),"\n",
    "LASSO Updated Mean Estimate of of phi_21 across all the simulations                               =",mean(exp(d$lasso.updated.theta[5,])),"\n",
    "LASSO Updated Mean Estimate of of phi_22 across all the simulations                               =",mean(exp(d$lasso.updated.theta[6,])),"\n",
    "LASSO Updated Mean Estimate of of phi_30 across all the simulations                               =",mean(exp(d$lasso.updated.theta[7,])),"\n",
    "LASSO Updated Mean Estimate of of phi_31 across all the simulations                               =",mean(exp(d$lasso.updated.theta[8,])),"\n",
    "LASSO Updated Mean Estimate of of phi_32 across all the simulations                               =",mean(exp(d$lasso.updated.theta[9,])),"\n",
    "LASSO Updated Mean Estimate of of phi_33 across all the simulations                               =",mean(exp(d$lasso.updated.theta[10,])),"\n",
    "LASSO Updated Mean Estimate of of frailty param across all the simulations                        =",mean(exp(d$lasso.updated.theta[11,])),"\n\n\n\n"
    
)

cat("LASSO Updated Standard deviation of of phi_10 across all the simulations                               =",sd(exp(d$lasso.updated.theta[1,])),"\n",
    "LASSO Updated Standard deviation of of phi_11 across all the simulations                               =",sd(exp(d$lasso.updated.theta[2,])),"\n",
    "LASSO Updated Standard deviation of of phi_12 across all the simulations                               =",sd(exp(d$lasso.updated.theta[3,])),"\n",
    "LASSO Updated Standard deviation of of phi_20 across all the simulations                               =",sd(exp(d$lasso.updated.theta[4,])),"\n",
    "LASSO Updated Standard deviation of of phi_21 across all the simulations                               =",sd(exp(d$lasso.updated.theta[5,])),"\n",
    "LASSO Updated Standard deviation of of phi_22 across all the simulations                               =",sd(exp(d$lasso.updated.theta[6,])),"\n",
    "LASSO Updated Standard deviation of of phi_30 across all the simulations                               =",sd(exp(d$lasso.updated.theta[7,])),"\n",
    "LASSO Updated Standard deviation of of phi_31 across all the simulations                               =",sd(exp(d$lasso.updated.theta[8,])),"\n",
    "LASSO Updated Standard deviation of of phi_32 across all the simulations                               =",sd(exp(d$lasso.updated.theta[9,])),"\n",
    "LASSO Updated Standard deviation of of phi_33 across all the simulations                               =",sd(exp(d$lasso.updated.theta[10,])),"\n",
    "LASSO Updated Standard deviation of of frailty param across all the simulations                        =",sd(exp(d$lasso.updated.theta[11,])),"\n\n\n\n"
    
)


cat("LASSO Updated Bias Estimate of theta (frailty parameter) across all simulations                 =",mean(exp(d$lasso.updated.theta[11,]))-theta.given,"\n\n\n"
)



#UPDATED ESTIMATES of  NUISANCE PARAMS BASED on **** ADAPTIVE LASSO:
#on LOG scale:
cat("Estimation results using ADAPTIVE LASSO Penalty based on the log scale:", "\n \n")

cat("Adaptive LASSO Updated Mean Estimate of log of phi_10 across all the simulations                               =",mean(d$adalasso.updated.theta[1,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_11 across all the simulations                               =",mean(d$adalasso.updated.theta[2,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_12 across all the simulations                               =",mean(d$adalasso.updated.theta[3,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_20 across all the simulations                               =",mean(d$adalasso.updated.theta[4,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_21 across all the simulations                               =",mean(d$adalasso.updated.theta[5,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_22 across all the simulations                               =",mean(d$adalasso.updated.theta[6,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_30 across all the simulations                               =",mean(d$adalasso.updated.theta[7,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_31 across all the simulations                               =",mean(d$adalasso.updated.theta[8,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_32 across all the simulations                               =",mean(d$adalasso.updated.theta[9,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of phi_33 across all the simulations                               =",mean(d$adalasso.updated.theta[10,]),"\n",
    "AdaptiveLASSO Updated Mean Estimate of log of frailty param across all the simulations                        =",mean(d$adalasso.updated.theta[11,]),"\n\n\n"
    
)


cat("AdaptiveLASSO Updated Standard deviation of log of phi_10 across all the simulations                               =",sd(d$adalasso.updated.theta[1,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_11 across all the simulations                               =",sd(d$adalasso.updated.theta[2,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_12 across all the simulations                               =",sd(d$adalasso.updated.theta[3,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_20 across all the simulations                               =",sd(d$adalasso.updated.theta[4,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_21 across all the simulations                               =",sd(d$adalasso.updated.theta[5,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_22 across all the simulations                               =",sd(d$adalasso.updated.theta[6,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_30 across all the simulations                               =",sd(d$adalasso.updated.theta[7,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_31 across all the simulations                               =",sd(d$adalasso.updated.theta[8,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_32 across all the simulations                               =",sd(d$adalasso.updated.theta[9,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of phi_33 across all the simulations                               =",sd(d$adalasso.updated.theta[10,]),"\n",
    "AdaptiveLASSO Updated Standard deviation of log of frailty param across all the simulations                        =",sd(d$adalasso.updated.theta[11,]),"\n\n\n\n"
    
)

cat("AdaptiveLASSO Updated Bias Estimate of log of theta (frailty parameter) across all the simulations             =",mean(d$adalasso.updated.theta[11,])-log(theta.given),"\n\n\n"
)







#on ORIGINAL scale:

cat("ADAPTIVELASSO Updated Mean Estimate of of phi_10 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[1,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_11 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[2,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_12 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[3,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_20 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[4,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_21 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[5,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_22 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[6,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_30 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[7,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_31 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[8,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_32 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[9,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of phi_33 across all the simulations                               =",mean(exp(d$adalasso.updated.theta[10,])),"\n",
    "ADAPTIVELASSO Updated Mean Estimate of of frailty param across all the simulations                        =",mean(exp(d$adalasso.updated.theta[11,])),"\n\n\n\n"
    
)

cat("ADAPTIVELASSO Updated Standard deviation of of phi_10 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[1,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_11 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[2,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_12 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[3,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_20 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[4,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_21 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[5,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_22 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[6,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_30 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[7,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_31 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[8,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_32 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[9,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of phi_33 across all the simulations                               =",sd(exp(d$adalasso.updated.theta[10,])),"\n",
    "ADAPTIVELASSO Updated Standard deviation of of frailty param across all the simulations                        =",sd(exp(d$adalasso.updated.theta[11,])),"\n\n\n\n"
    
)


cat("ADAPTIVELASSO Updated Bias Estimate of theta (frailty parameter) across all simulations                 =",mean(exp(d$adalasso.updated.theta[11,]))-theta.given,"\n\n\n"
)



cat("start time                                                                            =",starttime,"\n",
    "end time                                                                              =",endtime,"\n\n\n",
    "Simulation Result for rho                                                             =",rho2,"\n\n\n",
    "dimension after cutting left-trunc                                                    =",d$dimension.final.mean,"\n",
    "Mean of cens 00                                                                       =",d$censoring.rate.final[1],"\n",
    "Mean of obs 01                                                                        =",d$censoring.rate.final[2],"\n",
    "Mean of cens 10                                                                       =",d$censoring.rate.final[3],"\n",
    "Mean of obs 11                                                                        =",d$censoring.rate.final[4],"\n",
    "MSE.oracle                                                                            =",median(d$MSE.calculator.oracle[!is.na(d$MSE.calculator.oracle)]),"\n\n\n",
    "SD of MSE oracle                                                                      =",d$SD.oracle,"\n\n\n",
    
    "BAR:\n",
    "TP in bar                                                                             =",mean(d$TP.bar[!is.na(d$TP.bar)]),"\n",
    "FP in bar                                                                             =",mean(d$FP.bar[!is.na(d$FP.bar)]),"\n",
    "FN in bar                                                                             =",mean(d$FN.bar[!is.na(d$FN.bar)]),"\n",
    "MCV in bar                                                                            =",d$MC.final.bar,"\n",
    "MMSE in bar                                                                           =",median(d$MSE.bar[!is.na(d$MSE.bar)]),"\n",
    "SD of MSE in bar                                                                      =",d$SD.bar,"\n",
    "# of unsucessful iterations                                                           =",sum(is.na(d$TP.bar)),"\n\n\n",
    
    "LASSO:\n",
    "TP in lasso                                                                           =",mean(d$TP.lasso[!is.na(d$TP.lasso)]),"\n",
    "FP in lasso                                                                           =",mean(d$FP.lasso[!is.na(d$FP.lasso)]),"\n",
    "FN in lasso                                                                           =",mean(d$FN.lasso[!is.na(d$FN.lasso)]),"\n",
    "MCV in lasso                                                                          =",d$MC.final.lasso,"\n",
    "MMSE in lasso                                                                         =",median(d$MSE.lasso[!is.na(d$MSE.lasso)]),"\n",
    "SD of MSE in lasso                                                                    =",d$SD.lasso,"\n",
    "# of unsucessful iterations                                                           =",sum(is.na(d$TP.lasso)),"\n\n\n",
    
    "Adaptive lasso\n",
    "TP in .ada                                                                            =",mean(d$TP.ada[!is.na(d$TP.ada)]),"\n",
    "FP in .ada                                                                            =",mean(d$FP.ada[!is.na(d$FP.ada)]),"\n",
    "FN in ada                                                                             =",mean(d$FN.ada[!is.na(d$FN.ada)]),"\n",
    "MCV in .ada                                                                           =",mean(d$MC.final.ada),"\n",
    "MMSE in .ada                                                                          =",median(d$MSE.ada[!is.na(d$MSE.ada)]),"\n",
    "SD of MSE in .ada                                                                     =",d$SD.ada,"\n",
    "# of unsucessful iterations                                                           =",sum(is.na(d$TP.ada)),"\n\n\n\n\n\n\n\n\n")

cat("theta.bar                    =",d$theta.bar.mean,"\n",
    "theta.ada                    =",d$theta.ada.mean,"\n",
    "theta.lasso                  =",d$theta.lasso.mean,"\n",
    "SD of theta[7]               =",sd(d$theta.bar[7,]),"\n\n\n\n\n\n\n\n\n"
)

cat("complete vector of theta[1]          =",d$theta.bar[1,],"\n",
    "complete vector of theta[2]          =",d$theta.bar[2,],"\n",
    "complete vector of theta[3]          =",d$theta.lasso[3,],"\n",
    "complete vector of theta[4]          =",d$theta.lasso[4,],"\n",
    "complete vector of theta[5]          =",d$theta.lasso[5,],"\n",
    "complete vector of theta[6]          =",d$theta.lasso[6,],"\n",
    "complete vector of theta[7]          =",d$theta.lasso[7,],"\n",
    "complete vector of theta[8]          =",d$theta.lasso[8,],"\n",
    "complete vector of theta[9]          =",d$theta.lasso[9,],"\n",
    "complete vector of theta[10]          =",d$theta.lasso[10,],"\n",
    "complete vector of theta[11]          =",d$theta.lasso[11,],"\n\n\n\n\n\n\n\n\n\n")

cat("MSE's for rho                    =",rho2,"\n\n\n\n\n\n\n",
    "MSE for BAR                      =",d$MSE.bar,"\n\n",
    "MSE for Lasso                    =",d$MSE.lasso,"\n\n",
    "MSE for Adaptive Lasso           =",d$MSE.ada,"\n\n\n\n\n\n\n")
cat(" rho                                   =",rho2,"\n\n",
    "Mean MSE for BAR                      =",mean(d$MSE.bar),"\n\n",
    "Mean MSE for Lasso                    =",mean(d$MSE.lasso),"\n\n",
    "Mean MSE for Adaptive Lasso           =",mean(d$MSE.ada),"\n\n\n\n\n\n\n")

cat(" rho                                     =",rho2,"\n\n",
    "Median MSE for BAR                      =",median(d$MSE.bar),"\n\n",
    "Median MSE for Lasso                    =",median(d$MSE.lasso),"\n\n",
    "Median MSE for Adaptive Lasso           =",median(d$MSE.ada),"\n\n\n\n\n\n\n")

# print("mean bar estimates:")
# apply(d$beta.selected.bar,1,mean)
print("mean bar estimates:")
beta.selected.adjusted=d$beta.selected.bar
beta.selected.adjusted[abs(beta.selected.adjusted)<tol]=0
apply(beta.selected.adjusted,1,mean)
print("mean ada estimates:")
apply(d$beta.selected.ada,1,mean)
print("mean lasso estimates:")
apply(d$beta.selected.lasso,1,mean)

print("selection frequencies for BAR:")
rowSums(beta.selected.adjusted!=0)

print("selection frequencies for ADALASSO:")
rowSums(d$beta.selected.ada!=0)

print("selection frequencies for LASSO:")
rowSums(d$beta.selected.lasso!=0)



#UNPEN
cat("complete vectors because I need to run batch by batch so that we can get SD and Mean for different size of batches:","\n\n")
print("complete vector of phi_10 for unpen:")
cat(d$theta.bar[1,],"\n\n")
print("complete vector of phi_11 for unpen:")
cat(d$theta.bar[2,],"\n\n")
print("complete vector of phi_12 for unpen:")
cat(d$theta.bar[3,],"\n\n")
print("complete vector of phi_20 for unpen:")
cat(d$theta.bar[4,],"\n\n")
print("complete vector of phi_21 for unpen:")
cat(d$theta.bar[5,],"\n\n")
print("complete vector of phi_22 for unpen:")
cat(d$theta.bar[6,],"\n\n")
print("complete vector of phi_30 for unpen:")
cat(d$theta.bar[7,],"\n\n")

print("complete vector of phi_31 for unpen:")
cat(d$theta.bar[8,],"\n\n")

print("complete vector of phi_32 for unpen:")
cat(d$theta.bar[9,],"\n\n")

print("complete vector of phi_33 for unpen:")
cat(d$theta.bar[10,],"\n\n")

print("complete vector of frailty param for unpen:")
cat(d$theta.bar[11,],"\n\n\n\n\n\n")







#BAR:

print("complete vector of phi_10 for BAR:")
cat(d$BAR.updated.theta[1,],"\n\n")
print("complete vector of phi_11 for BAR:")
cat(d$BAR.updated.theta[2,],"\n\n")
print("complete vector of phi_12 for BAR:")
cat(d$BAR.updated.theta[3,],"\n\n")
print("complete vector of phi_20 for BAR:")
cat(d$BAR.updated.theta[4,],"\n\n")
print("complete vector of phi_21 for BAR:")
cat(d$BAR.updated.theta[5,],"\n\n")
print("complete vector of phi_22 for BAR:")
cat(d$BAR.updated.theta[6,],"\n\n")
print("complete vector of phi_30 for BAR:")
cat(d$BAR.updated.theta[7,],"\n\n")

print("complete vector of phi_31 for BAR:")
cat(d$BAR.updated.theta[8,],"\n\n")

print("complete vector of phi_32 for BAR:")
cat(d$BAR.updated.theta[9,],"\n\n")

print("complete vector of phi_33 for BAR:")
cat(d$BAR.updated.theta[10,],"\n\n")

print("complete vector of frailty param for BAR:")
cat(d$BAR.updated.theta[11,],"\n\n\n\n\n\n")



#LASSO:

print("complete vector of phi_10 for LASSO:")
cat(d$lasso.updated.theta[1,],"\n\n")
print("complete vector of phi_11 for LASSO:")
cat(d$lasso.updated.theta[2,],"\n\n")
print("complete vector of phi_12 for LASSO:")
cat(d$lasso.updated.theta[3,],"\n\n")
print("complete vector of phi_20 for LASSO:")
cat(d$lasso.updated.theta[4,],"\n\n")
print("complete vector of phi_21 for LASSO:")
cat(d$lasso.updated.theta[5,],"\n\n")
print("complete vector of phi_22 for LASSO:")
cat(d$lasso.updated.theta[6,],"\n\n")
print("complete vector of phi_30 for LASSO:")
cat(d$lasso.updated.theta[7,],"\n\n")

print("complete vector of phi_31 for LASSO:")
cat(d$lasso.updated.theta[8,],"\n\n")

print("complete vector of phi_32 for LASSO:")
cat(d$lasso.updated.theta[9,],"\n\n")

print("complete vector of phi_33 for LASSO:")
cat(d$lasso.updated.theta[10,],"\n\n")

print("complete vector of frailty param for LASSO:")
cat(d$lasso.updated.theta[11,],"\n\n\n\n\n\n")



#ADALASSO:

print("complete vector of phi_10 for ADALASSO:")
cat(d$adalasso.updated.theta[1,],"\n\n")
print("complete vector of phi_11 for ADALASSO:")
cat(d$adalasso.updated.theta[2,],"\n\n")
print("complete vector of phi_12 for ADALASSO:")
cat(d$adalasso.updated.theta[3,],"\n\n")
print("complete vector of phi_20 for ADALASSO:")
cat(d$adalasso.updated.theta[4,],"\n\n")
print("complete vector of phi_21 for ADALASSO:")
cat(d$adalasso.updated.theta[5,],"\n\n")
print("complete vector of phi_22 for ADALASSO:")
cat(d$adalasso.updated.theta[6,],"\n\n")
print("complete vector of phi_30 for ADALASSO:")
cat(d$adalasso.updated.theta[7,],"\n\n")

print("complete vector of phi_31 for ADALASSO:")
cat(d$adalasso.updated.theta[8,],"\n\n")

print("complete vector of phi_32 for ADALASSO:")
cat(d$adalasso.updated.theta[9,],"\n\n")

print("complete vector of phi_33 for ADALASSO:")
cat(d$adalasso.updated.theta[10,],"\n\n")

print("complete vector of frailty param for ADALASSO:")
cat(d$adalasso.updated.theta[11,],"\n\n\n\n\n\n")
print("this is d$MSE.calculator.oracle:")
print(d$MSE.calculator.oracle)
print("This is d$MSE.bar")
print(d$MSE.bar)
print("This is d$MSE.lasso")
print(d$MSE.lasso)
print("This is d$MSE.ada")
print(d$MSE.ada)
sink()



