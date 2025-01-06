###Revised Program for FusedSCAD using Reeder's new R package with "bs" hazard
###to select varibles in a semiparametric semi-competing risks model;
### The data simulation functions are revised in 
###"sim-functions-weibull-fatemehrevised-Dec10.R", gamma.true were added , 
###since Reeder's code needs gamma.true, but not sure how to apply the program to real data
### because gamma.true is not observed. This needs further investigation. 
### Noted on Dec 10, 2024.
### There are 3 cased for different q_n=12, 15, 16, one needs to manually change them
### We can run the program with 3 different frailty terms: Gamma (correct specified model),
### Log-Normal,  Log Gaussian Mixture model (misspecified models)
### Need include another R program for functions to generate simulated data;

cat("The program started: ", date(), "\n \n")
options(digits=5,width=200)
##R library needed to run this program
##Some commented lines can be removed late;

#library("Rcpp")
#library("RcppArmadillo")
#library("remotes")
#library("SemiCompRisks")
#library("SemiCompRisksPen") #This pacakage is replaced by library(SemiCompRisksFreq) 

# library(devtools)
# install_github(repo = "https://github.com/harrisonreeder/SemiCompRisksFreq")
library(SemiCompRisksFreq)  #Use this only, don't use  library(SemiCompRisksPen)
#available at https://github.com/harrisonreeder/SemiCompRisksFreq
library(knitr)

#run the following two lines to install semicompriskspen by H. Reeder:

# remotes::install_github("coatless-mac/macrtools")
# remotes::install_github("harrisonreeder/SemiCompRisksPen")
# library(SemiCompRisksPen)

#other ways I tried for installing the package:
#remotes::install_github("harrisonreeder/SemiCompRisks") #there is no package called "semicoprisks". It is "semicompriskspen" only, apparently. 
#devtools::install_local("/Users/fmahmoudi/Desktop/Reeder/SemiCompRisksExamples-main.zip")

randseed=2023
set.seed(randseed) 

##Change the following settings to get different simulation scenarios;

#Sample size:
n=300     #change beta too for different n; (c1, c2) for different censoring rate;
#function "genWB.simData.noseed.allcontin.(qn)cov" for different q_n;
#theta.given= for different frailty variance;
#frailty_distribution= for different failty distribution;


#Number of simulation iterations at B=500:
B=500

censoring="high"

nb=c(4,4,4)     #Number of basis functions in cubic BSplines,  
#4 is the minimum for 0 number of interior knots;
verb=0 #verbose: 0 or 1;
#Numeric indicating the amount of iteration information 


#To genearte One simulation, just write seeds=c(1). To have ten simulations, c(1:10), etc.
seeds=c(1)  #not used in this program.

#Here, we determine the frailty distribution. Could be gamma or LN. theta.given is the true parameter that we set. 
#theta.given is variance (sigma^2) for LN or variance for Gamma.








variance=0.25
#variance=1.00

# Log-scale GMM parameters for variance =0.25 and mean=1:
if (variance==0.25){
  mu1 <- -0.0735             # Mean of the first component on the log scale
  mu2 <- -0.1470            # Mean of the second component on the log scale
  sigma1 <- 0.3834  # Standard deviation of the first component on the log scale
  sigma2 <- 0.5422  # Standard deviation of the second component on the log scale
  
}else{
  mu1 <- -0.2228             # Mean of the first component on the log scale
  mu2 <- -0.4457           # Mean of the second component on the log scale
  sigma1 <- 0.6676 # Standard deviation of the first component on the log scale
  sigma2 <- 0.9441 # Standard deviation of the second component on the log scale
  
}

p1 <- 0.5            # Probability of sampling from the first component
p2 <- 0.5            # Probability of sampling from the second component

# Calculate the variance on the positive scale (just for verification)
variance_GM <- p1 * (exp(2 * mu1 + sigma1^2) * (exp(sigma1^2) - 1)) +
  p2 * (exp(2 * mu2 + sigma2^2) * (exp(sigma2^2) - 1))

# frailty_distribution="Gamma"
# frailty_distribution="LogNormal"
frailty_distribution="GaussianMixture"

#Here, we determine the frailty distribution. Could be gamma or LN. theta.given is the true parameter that we set. 
#theta.given is variance (sigma^2) for LN or variance for Gamma.
# theta.given=0.25 #We test two values, 0.25 and 1.
theta.given=variance_GM


#We test two values, 0.25 and 1.



theta.true=theta.given

#theta.true = ifelse(frailty == T, theta.given, 0) #From the other simu.. program;
#We need theta.true with phi.true and beta.true;
#theta.init= -1  #This actuallry gives an initial value for log(theta)
#Set a random theta.init around log(theta.given) in the "varsel-functions..-tes.R" code;
# e.g., theta.init=log(theta.given)+rnorm(1, mean=0, sd=0.1)

# frailty_distribution="Gamma"
# frailty_distribution="LogNormal"
frailty_distribution="GaussianMixture"




#?solution_path_function()
#this is the help function that has all the functions in it. Using this help and the example that Reeder has provided, you might 
#be able to find how to run his code and apply their mothed to your data.
##help(package="SemiCompRisksPen")
#??simID_PW

##getwd()
#setwd("/Users/xlu/Library/CloudStorage/OneDrive-UniversityofCalgary/XLu2020/Fatemeh-PhD-2019-2023/Fatemeh-PhD-Thesis-papers/Project2-SemicompetingRisks/SemiComp-Revision2-Aug3-Simu1")
##setwd("/Users/fmahmoudi/Desktop/Reeder/code/revised-mywork")
#source("sim-functions-weibull-fatemehrevised-april30.R")
# source("sim-functions-weibull-fatemehrevised-Dec10.R")  #run this for gamma and LN
source("simfuncforGM.R")  #run this for gaussian mixture
#XL added  gamma.true in output, but there were two frailty terms in the output, 
##need to change it, also need to change 
# source("varsel-functions-fatemeh-revised-apr30.R")
#functions for BAR simulation study, the fused methods don't need them;

########choose different beta-vectors for different sample size n #############
###################n=?,rho=?################

#See Table 2 in the manuscript, p_n should be q_n;
#true beta parameters:
#p_n=12: #When p_n=12, set n=100.



weibull.param.log=c(-4,0.18,-4,0.2,-11,1.7)

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
#Censoring rate=low:

if (censoring=="low"){
  c1=1000
  c2=1000
}
if (censoring=="high"){
  c1=45
  c2=45
}



frailty=TRUE

#Here, I have defined a way to have two types of beta parameters. Variables across three risks in the same locations or in different locations. 
#fullyshared is when they are in same location just like our work. 
#overlapping is for when they are in different palce. 
#Some changes required if working with overlapping though. Right now (August 2022), we focusing on the fullyshared only.
# risks.status="overlapping"
risks.status="fullyshared"

#grouping effect study or individual variable selection:
grp.or.indiv="indiv"
# grp.or.indiv="grp"

#Optimization method. Could b optim or nlm. Both work. nlm works better sometimes when working with bernstein polynomials.
method="optim"
# method="nlm"

yvalslength=500 # #of some y values that are required to be generated while approximating baseline hazards using bernstein polynomials. 
#parameters to change the left-trunctaion setting:
lt1=0
lt2=0.3

#Sample size:
#n=500

#Df of BP but this file is for Weibull. Do not chnage them here. There is no need to these:
m=c(2,2,3) #degree of Bernstein polynomial for three transitions corresponding to 3 transitions 1,2,and 3.

#myseed=paste(min(seeds),max(seeds),sep=":")

#rho is the correlation between Z_i and Z_j set to 0.5 in this work:
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

# not required:
# lambdavec.lasso=seq(8.9,9.4,0.05)
# lambdavec.ada=seq(8.9,9.4,0.05)
# lambdavec.bar=c(1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6)

#Tuning parameters: 
#Have found them experimentally:
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
optim_starting_values=matrix(0,length(seeds),7)
censoring.observ.rates=matrix(0,length(seeds),4)#first col: both censored, 2nd col:nonterminal not observed, 3rd col: terminal censored, 4th col: both observed.
MSE.calculator.oracle=vector()
p=length(beta1.true)
p.oracle=sum(beta1.true!=0)
dimension=vector()

unpen.est=matrix(NA,length(seeds),3*p)
unpen.est.oracle=matrix(NA,length(seeds),3*p.oracle)

censoring.rate=c()

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

logtheta.init=NULL

BAR.updated.theta=matrix(0,7,length(seeds))
TP.bar=vector()
FP.bar=vector()
FN.bar=vector()
MSE.bar=vector()
GCV.bar=vector()
theta.bar=matrix(0,length(lambdavec.bar),length(seeds))
GCV.selected.bar=matrix(0,length(lambdavec.bar),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.bar))
FP.for.that.lam=vector()
MSE.for.that.lam=vector()
theta.for.that.lam.bar=matrix(0,7,length(lambdavec.bar))
betavarsel.bar=matrix(0,3*p,length(lambdavec.bar))
beta.selected.bar=matrix(0,3*p,length(seeds))


adalasso.updated.theta=matrix(0,7,length(seeds))
TP.ada=vector()
FP.ada=vector()
FN.ada=vector()
MSE.ada=vector()
theta.ada=matrix(0,length(lambdavec.ada),length(seeds))
GCV.AdaLasso=vector()
GCV.selected.ada=matrix(0,length(lambdavec.ada),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.ada))
FP.for.that.lam=vector()
FN.for.that.lam=vector()
MSE.for.that.lam=vector()
theta.for.that.lam.ada=matrix(0,7,length(lambdavec.ada))
betavarsel.AdaLasso=matrix(0,3*p,length(lambdavec.ada))
beta.selected.ada=matrix(0,3*p,length(seeds))

lasso.updated.theta=matrix(0,7,length(seeds))
TP.lasso=vector()
FP.lasso=vector()
FN.lasso=vector()
MSE.lasso=vector()
theta.lasso=matrix(0,length(lambdavec.lasso),length(seeds))
GCV.lasso=vector()
GCV.selected.lasso=matrix(0,length(lambdavec.lasso),length(seeds))
TP.for.that.lam=rep(0,length(lambdavec.lasso))
FP.for.that.lam=vector()
FN.for.that.lam=vector()
MSE.for.that.lam=vector()
theta.for.that.lam.lasso=matrix(0,7,length(lambdavec.lasso))
betavarsel.lasso=matrix(0,3*p,length(lambdavec.lasso))
beta.selected.lasso=matrix(0,3*p,length(seeds))

### Ignore the above declaration,  use the following:
TP.nofused.vec=vector()
FP.nofused.vec=vector()
MCV.nofused.vec=vector()
MSE.nofused.vec=vector()
theta.nofused.vec=vector()

TP.fused.vec=vector()
FP.fused.vec=vector()
MCV.fused.vec=vector()
MSE.fused.vec=vector()
theta.fused.vec=vector()

TP.mle.vec=vector()
FP.mle.vec=vector()
MCV.mle.vec=vector()
MSE.mle.vec=vector()
theta.mle.vec=vector()

### Start simulation iteration here ######
for (j in 1:B){
  cat("This is the", j, "th iteration.", "\n")
  #generate a set of data for q_n=12: 
  #I am making it flexible to work for any n and pn and qn:
  
  if (n==100){
    WB.simData=genWB.simData.noseed.allcontin.12cov (mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n =n,rho=rho2, frailty = T,risks.status,theta.given,frailty_distribution)
    
  }
  if (n==300){
    WB.simData=genWB.simData.noseed.allcontin.15cov (mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n =n,rho=rho2, frailty = T,risks.status,theta.given,frailty_distribution)
    
  }
  if (n==500){
    WB.simData=genWB.simData.noseed.allcontin.16cov (mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n =n,rho=rho2, frailty = T,risks.status,theta.given,frailty_distribution)
    
  }
  #head(WB.simData)
  
  ## Use more function from "sim-functions-weibull-fatemehrevised-Dec10.R", 
  ## for different q_n=15, 16, ...., 
  ## Also need to change the parameter seting for beta in the begining.
  
  
  
  ###aplying Reeder code to my data set that was just generated:  -----------
  
  # library(SemiCompRisksPen)
  
  # p=(dim(WB.simData$Y)[2])-5
  p=(dim(WB.simData$Y)[2])-6    #We have included gamma.true
  
  Y=WB.simData$Y
  Ylength=dim(Y)[2]
  
  ##Y[,(6:Ylength)]->CovMat
  Y[,(7:Ylength)]->CovMat           #We have included gamma.true
  
  para=c(WB.simData$WBpara) 
  
  y1=WB.simData$Y$y1
  y2=WB.simData$Y$y2
  
  
  delta1=WB.simData$Y$delta1
  delta2=WB.simData$Y$delta2
  
  Xmat1=Xmat2=Xmat3=CovMat #dim: nx(q1+q2+q3) : nx(3*pn). 
  
  Xmat1=as.matrix(Xmat1)
  Xmat2=as.matrix(Xmat2)
  Xmat3=as.matrix(Xmat3)
  
  
  ###From here, XL deleted all the lines from the original program FusedLasso-TestAug15.R##
  ###Then copy and paste all the lines from Reeder's example program;
  
  ##Form simData like the following:
  #> head(simData)
  #       y1 delta1     y2 delta2 yL gamma.true       X1       X2       X3       X4       X5       X6        X7        # X8        X9      X10      X11      X12      X13      X14      X15       X16  X17      X18       X19       X20      X21       X22      X23      X24       X25
  
  yL=WB.simData$Y$L          # yL=rep(0, n) if no left-truncation, n is the sample size;
  gamma.true=WB.simData$Y$gamma.true  #XL added this in the output;
  
  Xmat=Xmat1
  nP=ncol(Xmat)
  colnames(Xmat) <- paste0("X",1:nP)                #change the column names;
  
  mysimData=data.frame(y1, delta1, y2, delta2, yL, gamma.true, Xmat)
  
  #cat("The program started: ", date(), "\n \n")
  
  # library(devtools)
  # install_github(repo = "https://github.com/harrisonreeder/SemiCompRisksFreq")
  # library(SemiCompRisksFreq) 
  #available at https://github.com/harrisonreeder/SemiCompRisksFreq
  
  
  ##ANALYSIS VARIABLES##
  ##******************##
  
  ##modeling parameters
  
  #SELECT ONE OF THESE FOR YOUR MODEL FIT!
  # hazard_temp <- "pw" # piecewise model fit
  # hazard_temp <- "wb" # weibull model fit
  hazard_temp <- "wb" # B-spline model fit
  # Aliases for these are "wb", "pw", "bs", and "rp" respectively.
  
  
  model_temp <- "semi-Markov"
  frailty_temp <- TRUE
  
  penalty_temp <- "scad"
  a_temp <- 3.7
  
  #vector of fusion lasso lambda values to examine at each grid point
  lambda_fusedcoef_vec <- c(0,0.015,0.03)
  
  # N_path_steps = 100, there will be 100 \lambda_1 values, then combined with the above 
  # 3 \lambada_2 values, the total path grid point is 300;
  
  n.path.steps = 100
  
  #vector of lambda values for main penalty can also be manually specified (decreasing magnitude),
  #but it is easier to have it generated automatically in the function below
  
  fuse_type <- "h2h3" #other option is "allfuse" but we'll focus on fusing just h2 and h3
  temp_lambda_fusedcoef <- if(fuse_type == "h2h3") cbind(0,0,lambda_fusedcoef_vec) else cbind(lambda_fusedcoef_vec,lambda_fusedcoef_vec,lambda_fusedcoef_vec)
  
  
  #'path' for fusion lasso smoothing method, this generally doesn't need to be changed...
  mu_smooth_fused_path <- c(1e-2,1e-3,1e-4,1e-5,1e-6)
  
  #  The Nesterov smoothing parameter μ 
  #   we first iterate proximal gradient descent to convergence with μ large (e.g., 10−2 in our
  #   applications), and then decrease μ and further iterate to convergence, and so on until μ is #  sufficiently small. 
  
  #max number of iterations of proximal gradient descent per lambda grid point
  maxit_temp <- 1200
  #define fusion tolerance (how close estimates must be to be considered 'fused' from a df perspective)
  fusiontol_temp <- 1e-3
  #define selection tolerance (how large estimates must be to be considered 'nonzero')
  selecttol_temp <- 1e-4
  
  #make a formula for the 'full' model
  form_string <- paste0(colnames(Xmat),collapse = "+")
  form_temp <- Formula::as.Formula(paste0("y1 + delta1 | y2 + delta2 ~ ", form_string," | ",form_string," | ",form_string))
  
  #now, run the penalized algorithm
  temp_fit <- FreqID_HReg_Rpath(
    Formula=form_temp, data = mysimData, na.action="na.fail", subset=NULL, weights=NULL,
    hazard=hazard_temp, #if hazard is set to weibull, knots_list and nP0 are ignored.
    #if you want to fit the exact simulated model, set hazard to piecewise and
    #set knots_list equal to the true knots. If you specify knots_list, nP0 is ignored.
    
    #knots_list=knots_list_true, #knots_list = NULL,
    
    #if you want to fit a model with knots selected based on the data
    #set knots_list = NULL and set nP0 the number of baseline parameters per transition
    #nP0 = c(4,4,4),  
    #*Actually, it is  temp_fit$nP0 [1] 4 4 3
    
    knots_list = NULL,
    # nP0 = c(4,4,3)+2,       #(6, 6, 5) basis functions 
    nP0 = nb,    #c(2,2,2)+2,         #(4, 4, 4) basis functions 
    # If less than 4, get Error in get_default_knots(y = if (tolower(model) == "semi-markov") (y2 -  : 
    #  Cubic B-Spline Specification must have at least 4 parameters in each hazard.
    
    ##XL: my experments shows if we use hazard="bs", then the knots_list becomes the 
    ## the interior knots, the total number of baseline parameters becomes nP0+4-2;
    ## If I set knots_list = NULL, nP0=c(4,4,3),  the  total number of baseline parameters 
    ## is still nP0;
    model=model_temp, frailty=frailty_temp,
    penalty=penalty_temp, a=a_temp, select_tol=selecttol_temp,
    #now, instead of manually specifying a lambda path, you can just specify
    #what the smallest/final lambda you're interested in is (lambda_target), and
    #the total number of steps on the solution path. Function automatically
    #chooses a reasonable starting lambda and defines the rest of the path from there.
    lambda_path=NULL, lambda_target=0.001, N_path_steps = n.path.steps,
    penalty_fusedcoef="fusedlasso", lambda_fusedcoef_path=temp_lambda_fusedcoef,
    mu_smooth_path=mu_smooth_fused_path, fusion_tol=fusiontol_temp,
    startVals=NULL, standardize = TRUE,
    fit_method="prox_grad", maxit=maxit_temp, extra_starts=0,
    conv_crit = "nll_pen_change", conv_tol=1e-6, verbose=verb)
  
  #verbose: 0 or 1;
  #Numeric indicating the amount of iteration information 
  #should be printed to the user. Higher #numbers provide more detailed information to user, #but will slow down the algorithm.
  
  #fitting unpenalized model takes many of the same arguments
  temp_mle <- FreqID_HReg2(Formula=form_temp, data = mysimData, na.action="na.fail", subset=NULL, weights=NULL,
                           hazard=hazard_temp,
                           #if you want to fit the exact simulated model, set knots_list equal to the true knots
                           # knots_list=knots_list_true, #knots_list = NULL,
                           #if you want to fit a model with knots selected based on the data
                           #set knots_list = NULL and specify the number of baseline parameters per transition
                           #nP0 = c(4,4,4),
                           
                           #  knots_list = NULL, 
                           #  nP0 = c(4,4,3)+2,      #(6, 6, 5) basis functions 
                           nP0 =nb,  #c(2,2,2)+2,         #(4, 4, 4) basis functions 
                           # If less than 4, get Error in get_default_knots(y = if (tolower(model) == "semi-markov") (y2 -  : 
                           #  Cubic B-Spline Specification must have at least 4 parameters in each hazard.
                           model=model_temp, frailty=frailty_temp)
  
  # Function FreqID_HReg2() is different in the two R packages: 
  # ??  FreqID_HReg2 search online Help pages:
  # SemiCompRisksFreq uses nP0= argument;
  # SemiCompRisksPen uses   p0_vec= argument;
  # Both Fit Parametric Frailty Illness-Death Model for Semi-Competing Risks Data
  
  # Use one *Freq package correctly, otherwise, Had this: 
  # Error in FreqID_HReg2(Formula = form_temp, data = mysimData, 
  # na.action = "na.fail",  : unused arguments (weights = NULL, nP0 = c(2, 2, 2) + 2)
  
  #########I deleted the plots lines#############
  #plot the 'solution path' for each transition
  #library(dplyr)
  #library(tidyr)
  #library(ggplot2)
  
  #compare the fitted parameters to the true parameters
  
  ic_frame <- cbind(temp_fit$info, temp_fit$ics) %>%
    mutate(lambda_fusedcoef3_fct = as.factor(lambda_fusedcoef3))
  ic_frame_nofused <- ic_frame %>% filter(lambda_fusedcoef3==0)
  
  lambda1_min <- ic_frame$lambda1[which.min(ic_frame$BIC_unique)]
  lambda_fusedcoef3_min <- ic_frame$lambda_fusedcoef3[which.min(ic_frame$BIC_unique)]
  lambda1_nofused_min <- ic_frame_nofused$lambda1[which.min(ic_frame_nofused$BIC_unique)]
  
  
  para_best_nofused_BIC <- as.matrix(temp_fit$ests)[temp_fit$info$lambda1==lambda1_nofused_min &
                                                      temp_fit$info$lambda_fusedcoef3==0,]
  para_best_BIC <- as.matrix(temp_fit$ests)[temp_fit$info$lambda1==lambda1_min &
                                              temp_fit$info$lambda_fusedcoef3==lambda_fusedcoef3_min,]
  
  # (only works if you fit piecewise constant model
  # with the same knots from the data generation in the fitted model,
  # otherwise parameters don't line up)
  
  ##XL: Yes,When I used hazard=”bs”, nP0 = c(4,4,4) actually are the number of basis functions,
  ##If I used knots_list, then the nunber of interior knots=nP0-2 for B-splines, then when degree=3 (cubic spline), the order=3+1=4, the total number of B-spline basis functions is Number of interior knots + order =(4,4,3)-2+4=(6,6,5), this is what I obtained from Reeder’s sample code. 
  ## Since the minum number of knots is 2 (i.e., zero interior knots), and the cubic spline has order of 4, then the minimum value for the number of B-spline basis functions is nP0=(4,4,4) if knots_list is not used;
  
  
  #cbind(true=para_true_PW,
  #      est_nofused=para_best_nofused_BIC,
  #      est_fused=para_best_BIC,
  #      est_mle=temp_mle$estimate)
  
  #For bs or any other methods including PW, we revised the output as follows:
  #options(digits=5, width=200)
  
  # cat("Summary of the estimation for all the parameters:", "\n \n")
  
  #true parameters used in generating the data:
  theta.true = ifelse(frailty == T, theta.given, 0) 
  #From the other simu.. program;
  #We need theta.true with phi.true and beta.true;
  #Example: p_n=12: #When p_n=12, set n=100.
  #beta1.true = c(-0.8,1,1,0.9,rep(0,8))
  #beta2.true = c(1,1,1,0.9,rep(0,8))
  #beta3.true = c(-1,1,0.9,1,rep(0,8))
  #weibull.param.log=c(-4,0.18,-4,0.2,-11,1.7)
  
  ## We define our true parameters used in generating the data, 
  ## see the beginning of the program
  
  weibull.param=exp(weibull.param.log)
  phi1_true=weibull.param[1:2]
  phi2_true=weibull.param[3:4]
  phi3_true=weibull.param[5:6]
  theta_true=theta.true
  beta1_true=beta1.true
  beta2_true=beta2.true
  beta3_true=beta3.true
  
  ## From Reeder: 
  ## para_true_PW <- c(phi1_true,phi2_true,phi3_true,log(theta_true), beta1_true,beta2_true, beta3_true)
  
  para_true_PW <- c(phi1_true,phi2_true,phi3_true, log(theta_true), beta1_true,beta2_true,beta3_true)
  
  #norig.phi<-length(c(knots_list_true[[1]], knots_list_true[[2]], knots_list_true[[3]]))
  #number of the original phi parameters set in the simulation data
  ## For semiparametric model, we don't have true phi knots except the distribution parameters;
  
  ## Need make a loop to repeat simlation for B=500 times;
  #  B=1
  # cat("Number of replication in simu study.  =", B, "\n")
  # cat("Sample size                                          =", n, "\n")
  # cat("Number of preditors in each transition=", nP, "\n")
  # cat("Total Number of preditors                    =", nP*3, "\n")
  # cat("Variance of the frailty term                   =", theta.true, "\n \n")
  
  norig.phi=length(c(phi1_true, phi2_true, phi3_true))
  nest.phi<- sum(temp_fit$nP0)    
  # number of estimated phi parameters when knots_list_true is used, it is  norig.phi+2, 
  # if knots_list_true is not used, then it is nP0;
  
  true.phi=para_true_PW[1:norig.phi]
  # cat("true phi:", "\n")
  # print(true.phi)
  
  phi.result<-data.frame(nofused.BIC.phi=para_best_nofused_BIC[1:nest.phi], fused.BIC.phi=para_best_BIC[1:nest.phi], mle.phi=temp_mle$estimate[1:nest.phi])
  # cat("etimated phi only:", "\n")
  # print(phi.result)
  
  true.ltheta=para_true_PW[norig.phi+1]
  
  nofused.BIC.ltheta=para_best_nofused_BIC[nest.phi+1]
  fused.BIC.ltheta=para_best_BIC[nest.phi+1]
  mle.ltheta=temp_mle$estimate[nest.phi+1]
  
  #ltheta.result<-data.frame(true.ltheta, nofused.BIC.ltheta, fused.BIC.ltheta, mle.ltheta)
  #cat("compared estimated ltheta with true ltheta (log(theta)):", "\n")
  #print(ltheta.result)
  
  #cat("compared estimated theta with true theta:", "\n")
  #print(exp(ltheta.result))
  
  #beta.result<-data.frame(true.beta,  nofused.BIC.beta, fused.BIC.beta, mle.beta)
  #cat("compared estimated beta with true beta:", "\n")
  #print(beta.result)
  
  
  E_ZZprime=(1/n)*(t(cbind(Xmat1,Xmat2,Xmat3))%*%cbind(Xmat1,Xmat2,Xmat3))
  ##Xmat has been centered;
  true.beta=para_true_PW[-(1:(norig.phi+1))]
  nofused.BIC.beta=para_best_nofused_BIC[-(1:(nest.phi+1))]
  fused.BIC.beta=para_best_BIC[-(1:(nest.phi+1))]
  mle.beta=temp_mle$estimate[-(1:(nest.phi+1))]
  
  true.P=length(which(true.beta!=0))      #Number of true positive
  nonzero.truth.index=which(true.beta!=0)
  zero.truth.index=which(true.beta==0)
  
  TP.nofused=length(which(abs(nofused.BIC.beta[nonzero.truth.index])>tol))     #0.0001
  FP.nofused=length(which(abs(nofused.BIC.beta[zero.truth.index])>tol))
  MCV.nofused=(true.P-TP.nofused)+ FP.nofused     #FN+FP
  MSE.nofused=t(as.matrix(nofused.BIC.beta-true.beta))%*%(E_ZZprime)%*%as.matrix(nofused.BIC.beta-true.beta)
  
  TP.fused=length(which(abs(fused.BIC.beta[nonzero.truth.index])>tol))     #0.0001
  FP.fused=length(which(abs(fused.BIC.beta[zero.truth.index])>tol))
  MCV.fused=(true.P-TP.fused)+ FP.fused     #FN+FP
  MSE.fused=t(as.matrix(fused.BIC.beta-true.beta))%*%(E_ZZprime)%*%as.matrix(fused.BIC.beta-true.beta)
  
  TP.mle=length(which(abs(mle.beta[nonzero.truth.index])>tol))     #0.0001
  FP.mle=length(which(abs(mle.beta[zero.truth.index])>tol))       #mle does not do selection, 
  MCV.mle=(true.P-TP.mle)+ FP.mle    #FN+FP                           #just report them, not use it. 
  MSE.mle=t(as.matrix(mle.beta-true.beta))%*%(E_ZZprime)%*%as.matrix(mle.beta-true.beta)
  
  
  TP.nofused.vec[j]= TP.nofused
  FP.nofused.vec[j]= FP.nofused
  MCV.nofused.vec[j]=MCV.nofused
  MSE.nofused.vec[j]=MSE.nofused
  theta.nofused.vec[j]=exp(nofused.BIC.ltheta)
  
  TP.fused.vec[j]= TP.fused
  FP.fused.vec[j]= FP.fused
  MCV.fused.vec[j]=MCV.fused
  MSE.fused.vec[j]=MSE.fused
  theta.fused.vec[j]=exp(fused.BIC.ltheta)
  
  TP.mle.vec[j]= TP.mle
  FP.mle.vec[j]= FP.mle
  MCV.mle.vec[j]=MCV.mle
  MSE.mle.vec[j]=MSE.mle
  theta.mle.vec[j]=exp(mle.ltheta)
  
} #end of j loop for simulation replications;
## Need make a loop to repeat simlation for B=500 times;

cat("Random seed for the simulation          =", randseed, "\n")
cat("frailty_distribution                               =", frailty_distribution, "\n")
cat("Variance of the frailty term                  =", theta.true, "\n ")

cat("Censoring parameters (c1, c2)             =", c(c1,c2), "\n ")
cat("Tolerance value for sparsity                 =", tol, "\n ")
#grouping effect study or individual variable selection: grp.or.indiv="indiv"
#grouping effect study or individual variable selection:
cat("The study is group or individual selection study =", grp.or.indiv, "\n ")
cat("Correlaton parameter rho in rho^(abs(i-j))    =", rho2, "\n")

cat("Number of basis functions in cubic BSplines=" , nb, "\n")
cat("Number of replication in simu study   =", B, "\n")
cat("Sample size                                          =", n, "\n")
cat("Number of preditors in each transition=", nP, "\n")
cat("Total Number of preditors                    =", nP*3, "\n \n")


mean.TP.nofused=mean(TP.nofused.vec)
mean.FP.nofused=mean(FP.nofused.vec)
mean.MCV.nofused=mean(MCV.nofused.vec)
med.MSE.nofused=median(MSE.nofused.vec)
SD.MSE.nofused=sd(MSE.nofused.vec)
mean.theta.nofused=mean(theta.nofused.vec)
SD.theta.nofused=sd(theta.nofused.vec)

mean.TP.fused=mean(TP.fused.vec)
mean.FP.fused=mean(FP.fused.vec)
mean.MCV.fused=mean(MCV.fused.vec)
med.MSE.fused=median(MSE.fused.vec)
SD.MSE.fused=sd(MSE.fused.vec)
mean.theta.fused=mean(theta.fused.vec)
SD.theta.fused=sd(theta.fused.vec)

mean.TP.mle=mean(TP.mle.vec)
mean.FP.mle=mean(FP.mle.vec)
mean.MCV.mle=mean(MCV.mle.vec)
med.MSE.mle=median(MSE.mle.vec)
SD.MSE.mle=sd(MSE.mle.vec)
mean.theta.mle=mean(theta.mle.vec)
SD.theta.mle=sd(theta.mle.vec)


Result.nofused=c(mean.TP.nofused,mean.FP.nofused,mean.MCV.nofused,med.MSE.nofused,SD.MSE.nofused,mean.theta.nofused,SD.theta.nofused)
Result.fused=c(mean.TP.fused,mean.FP.fused,mean.MCV.fused,med.MSE.fused,SD.MSE.fused,mean.theta.fused,SD.theta.fused)
Result.mle=c(mean.TP.mle,mean.FP.mle,mean.MCV.mle,med.MSE.mle,SD.MSE.mle,mean.theta.mle,SD.theta.mle)

beta.true=data.frame(beta1_true,beta2_true,beta3_true)
phi.true=data.frame(phi1_true,phi2_true,phi3_true)

Result.all=rbind(Result.nofused, Result.fused, Result.mle)
colnames(Result.all) <- c("TP", "FP", "MCV", "MMSE", "SD(MMSE)", "hat.theta", "SD(hat.theta)") 
#paste0("X",1:nP)                #change the column names;
rownames(Result.all)<-c("Nofused", "Fused", "MLE")

cat("Three true phi vectors :",  "\n")
print(phi.true)
cat("Three true beta vectors:",  "\n")
print(beta.true)

cat("Result for nofused, fused and mle by Reeder's method:", "\n")
#print(Result.all)
print(knitr::kable(Result.all))

cat("The program ended:", date(), "\n \n")





