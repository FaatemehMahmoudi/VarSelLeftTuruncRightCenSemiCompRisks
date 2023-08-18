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
source("sim-functions-weibull.R")
source("varsel-functions-weibull.R")
###########################################################################
###################n=?,rho=?################
# options(digits=3,width=200)
#10 variables considered for each transition with the following values:
seeds=c(1:100)

beta1.true = c(-0.8,1,1,0.9,rep(0,8))
beta2.true = c(1,1,1,0.9,rep(0,8))
beta3.true = c(-1,1,0.9,1,rep(0,8))
weibull.param.log=c(-4,0.18,-4,0.2,-11,1.7,-1.4)
#Censoring rate=low:c1=c2=1000/high:c1=c2=45 corresponding to ~50%, and 65% of right-censored cases.
c1=1000
c2=1000
weibull.param.log=c(-4,0.18,-4,0.2,-11,1.7,-1.4)
#Lee et al parameters originally used in their paper:
# weibull.param.log=c(-9.98,1.05,-10.01,1.15,-5.92,0.92)
# c1=30
# c2=30
#65-65:
# lt1=0
#78-65:
# lt2=13
frailty=TRUE
# risks.status="overlapping"
risks.status="fullyshared"
#grouping effect study or individual variable selection:
grp.or.indiv="indiv"

method="optim"
yvalslength=500 # #of some y values that are required to be generated while approximating baseline hazards using bernstein polynomials. 
lt1=0
lt2=0.3
n=100
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

#indiv:
lambdavec.lasso=seq(8,9,0.15)
lambdavec.ada=seq(8,9,0.15)
# lambdavec.lasso=seq(5,6.5,1)
# lambdavec.ada=seq(2,3.5,1)
# lambdavec.bar=seq(1,2,0.25)
# lambdavec.bar=seq(1,2.5,0.25)
lambdavec.bar=seq(0.5,1.5,0.15)

#grp:
# lambdavec.lasso=seq(5,6,0.25)
# lambdavec.ada=seq(4,5,0.25)
# lambdavec.bar=seq(0.75,1.75,0.25)

xi=10
tol=1e-04#criterion for selecting zero and nonzero in BAR.

censoring.observ.rates=matrix(0,length(seeds),4)#first col: both censored, 2nd col:nonterminal not observed, 3rd col: terminal censored, 4th col: both observed.
MSE.calculator.oracle=vector()
p=length(beta1.true)
p.oracle=sum(beta1.true!=0)
dimension=vector()

unpen.est=matrix(NA,length(seeds),3*p)
unpen.est.oracle=matrix(NA,length(seeds),3*p.oracle)

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

TP.bar=vector()
FP.bar=vector()
FN.bar=vector()
MSE.bar=vector()
GCV.bar=vector()
GCV.selected.bar=matrix(0,length(lambdavec.bar),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.bar))
FP.for.that.lam=vector()
MSE.for.that.lam=vector()
betavarsel.bar=matrix(0,3*p,length(lambdavec.bar))
beta.selected.bar=matrix(0,3*p,length(seeds))

TP.ada=vector()
FP.ada=vector()
FN.ada=vector()
MSE.ada=vector()
GCV.AdaLasso=vector()
GCV.selected.ada=matrix(0,length(lambdavec.ada),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.ada))
FP.for.that.lam=vector()
FN.for.that.lam=vector()
MSE.for.that.lam=vector()
betavarsel.AdaLasso=matrix(0,3*p,length(lambdavec.ada))
beta.selected.ada=matrix(0,3*p,length(seeds))

TP.lasso=vector()
FP.lasso=vector()
FN.lasso=vector()
MSE.lasso=vector()
GCV.lasso=vector()
GCV.selected.lasso=matrix(0,length(lambdavec.lasso),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.lasso))
FP.for.that.lam=vector()
FN.for.that.lam=vector()
MSE.for.that.lam=vector()
betavarsel.lasso=matrix(0,3*p,length(lambdavec.lasso))
beta.selected.lasso=matrix(0,3*p,length(seeds))

starttime=date()
d=repeat.get.last.result(seeds,tol,n,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,weibull.param.log,lambdavec.lasso,lambdavec.ada,lambdavec.bar,xi,rho2,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli)

# d=repeat.get.last.result.bp(seeds,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho2,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,lambdavec.bar,lambdavec.lasso,lambdavec.ada)
endtime=date()

filename=filename="semi-p10-n300-c20-sg-r05"

sink(filename)

cat("n             =",n,"\n",
    "censoring     =",c1,"\n",
    "pn            =",length(beta1.true),"\n",
    "grporindiv    =",grp.or.indiv,"\n",
    "strongorweak  =",sum(beta1.true),"\n\n\n\n\n")

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

cat("MSE's for rho                    =",rho2,"\n\n\n\n\n\n\n",
    "MSE for BAR                      =",d$MSE.bar,"\n\n",
    "MSE for Lasso                    =",d$MSE.lasso,"\n\n",
    "MSE for Adaptive Lasso           =",d$MSE.ada,"\n\n\n\n\n\n\n")

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


print("sum of estimated beta coefficients row-wise for BAR:")
rowSums(beta.selected.adjusted)
print("sum of estimated beta coefficients row-wise for ADALASSO:")
rowSums(beta.selected.ada)
print("sum of estimated beta coefficients row-wise for LASSO:")
rowSums(beta.selected.lasso)


sink()



