library(dplyr)
library(survival)
require(stats)
# library(splines2)
library(pracma)  
# library(caret)
library(MASS)
# library(glmnet)
# library(lassoshooting)
WB.haz <- function(log.kappa, log.alpha, t){
  kappa = exp(log.kappa); alpha = exp(log.alpha)
  h = alpha * kappa * t^(alpha - 1)
  return(h)
}


simID.LT.noseed <- function(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                            alpha1.true, alpha2.true, alpha3.true,
                            kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt.type = "unif", lt,risks.status)
{
  
  n <- dim(x1)[1]
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]
  p3 <- dim(x3)[2]
  
  if(theta.true >0)
  {
    gamma.true <- rgamma(n, 1/theta.true, 1/theta.true)
  }
  if(theta.true == 0)
  {
    gamma.true <- rep(1, n)
  }
  
  LP1	<- as.vector(beta1.true %*% t(x1))
  LP2	<- as.vector(beta2.true %*% t(x2))
  LP3	<- as.vector(beta3.true %*% t(x3))
  
  Rind <- NULL
  R <- rweibull(n, shape = alpha1.true, scale = exp(-(log(kappa1.true) +
                                                        LP1 + log(gamma.true))/alpha1.true))
  D <- rweibull(n, shape = alpha2.true, scale = exp(-(log(kappa2.true) +
                                                        LP2 + log(gamma.true))/alpha2.true))
  
  yesR <- R < D
  
  D[yesR] <- R[yesR] + rweibull(sum(yesR), shape = alpha3.true,
                                scale = exp(-(log(kappa3.true) + LP3[yesR] + log(gamma.true[yesR]))/alpha3.true))
  delta1 <- rep(NA, n)
  delta2 <- rep(NA, n)
  y1 <- R
  y2 <- D
  Cen <- runif(n, cens[1], cens[2])
  ind01 <- which(D < R & D < Cen)
  y1[ind01] <- D[ind01]
  delta1[ind01] <- 0
  delta2[ind01] <- 1
  ind10 <- which(R < D & R < Cen & D >= Cen)
  y2[ind10] <- Cen[ind10]
  delta1[ind10] <- 1
  delta2[ind10] <- 0
  ind00 <- which(R >= Cen & D >= Cen)
  y1[ind00] <- Cen[ind00]
  y2[ind00] <- Cen[ind00]
  delta1[ind00] <- 0
  delta2[ind00] <- 0
  ind11 <- which(R < Cen & D < Cen & R < D)
  delta1[ind11] <- 1
  delta2[ind11] <- 1
  
  if (lt.type == "unif") L <- runif(n, lt[1], lt[2])
  if (lt.type == "norm") L <- rnorm(n, lt[1], lt[2])
  ret <- data.frame(cbind(y1, delta1, y2, delta2, L))
  return(ret)
}


genWB.simData.noseed.allcontin<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,10)
  for(k in 1:n){
    Sigma<-matrix(rep(0,10*10),10,10)
    for(i in 11:20){
      for(j in 11:20){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 10), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:15){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  
genWB.simData.noseed.allcontin.20cov<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth+twentyth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth+twentyth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth+twentyth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,20)
  for(k in 1:n){
    Sigma<-matrix(rep(0,20*20),20,20)
    for(i in 11:30){
      for(j in 11:30){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 20), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  thirteenth=as.vector(cov.mat[,13])
  fourteenth=as.vector(cov.mat[,14])
  fifteenth=as.vector(cov.mat[,15])
  sixteenth=as.vector(cov.mat[,16])
  seventeenth=as.vector(cov.mat[,17])
  eighteenth=as.vector(cov.mat[,18])
  ninteenth=as.vector(cov.mat[,19])
  twentyth=as.vector(cov.mat[,20])
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth,fifteenth,sixteenth,seventeenth,eighteenth,ninteenth,twentyth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:15){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  
genWB.simData.noseed.allcontin.12cov<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,12)
  for(k in 1:n){
    Sigma<-matrix(rep(0,12*12),12,12)
    for(i in 11:22){
      for(j in 11:22){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 12), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  
  data0=data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:17){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  
genWB.simData.noseed.allcontin.14cov<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,14)
  for(k in 1:n){
    Sigma<-matrix(rep(0,14*14),14,14)
    for(i in 11:24){
      for(j in 11:24){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 14), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  thirteenth=as.vector(cov.mat[,13])
  fourteenth=as.vector(cov.mat[,14])
  
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:19){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  

genWB.simData.noseed.allcontin.15cov<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,15)
  for(k in 1:n){
    Sigma<-matrix(rep(0,15*15),15,15)
    for(i in 11:25){
      for(j in 11:25){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 15), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  thirteenth=as.vector(cov.mat[,13])
  fourteenth=as.vector(cov.mat[,14])
  fifteenth=as.vector(cov.mat[,15])
  # sixteenth=as.vector(cov.mat[,16])
  # seventeenth=as.vector(cov.mat[,17])
  # eighteenth=as.vector(cov.mat[,18])
  # ninteenth=as.vector(cov.mat[,19])
  # twentyth=as.vector(cov.mat[,20])
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth,fifteenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:20){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  
genWB.simData.noseed.allcontin.18cov<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,18)
  for(k in 1:n){
    Sigma<-matrix(rep(0,18*18),18,18)
    for(i in 11:28){
      for(j in 11:28){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 18), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  thirteenth=as.vector(cov.mat[,13])
  fourteenth=as.vector(cov.mat[,14])
  fifteenth=as.vector(cov.mat[,15])
  sixteenth=as.vector(cov.mat[,16])
  seventeenth=as.vector(cov.mat[,17])
  eighteenth=as.vector(cov.mat[,18])
  
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth,fifteenth,sixteenth,seventeenth,eighteenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:23){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  
genWB.simData.noseed.allcontin.16cov<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,16)
  for(k in 1:n){
    Sigma<-matrix(rep(0,16*16),16,16)
    for(i in 11:26){
      for(j in 11:26){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 16), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  thirteenth=as.vector(cov.mat[,13])
  fourteenth=as.vector(cov.mat[,14])
  fifteenth=as.vector(cov.mat[,15])
  sixteenth=as.vector(cov.mat[,16])
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth,fifteenth,sixteenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:21){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  

genWB.simData.noseed.allcontin.19cov<- function(weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n , rho,frailty,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  library(MASS)
  cov.mat<-matrix(0,n,19)
  for(k in 1:n){
    Sigma<-matrix(rep(0,19*19),19,19)
    for(i in 11:29){
      for(j in 11:29){
        Sigma[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat[k,]<-mvrnorm(n = 1, rep(0, 19), Sigma)
  }
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  thirteenth=as.vector(cov.mat[,13])
  fourteenth=as.vector(cov.mat[,14])
  fifteenth=as.vector(cov.mat[,15])
  sixteenth=as.vector(cov.mat[,16])
  seventeenth=as.vector(cov.mat[,17])
  eighteenth=as.vector(cov.mat[,18])
  ninteenth=as.vector(cov.mat[,19])
  
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth,fifteenth,sixteenth,seventeenth,eighteenth,ninteenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:24){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  

genWB.simData.noseed.allcontin.grp<- function(G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n = n, rho=rho,frailty = T,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  cov.mat1<-matrix(0,n,length(G1))
  for(k in 1:n){
    Sigma1<-matrix(rep(0,length(G1)*length(G1)),length(G1),length(G1))
    for(i in 11:(10+length(G1))){
      for(j in 11:(10+length(G1))){
        Sigma1[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat1[k,]<-mvrnorm(n = 1, rep(0,length(G1)), Sigma1)
  }
  if (G1.distr=="normal"){
    cov.mat1=cov.mat1
  }else{
    if(G1.distr=="Bernoulli"){
      pvars1 <- pnorm(cov.mat1)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat1 <- qbinom(pvars1, 1, mu.bernoulli)}}
  #############second group:###########
  cov.mat2<-matrix(0,n,length(G2))
  for(k in 1:n){
    Sigma2<-matrix(rep(0,length(G2)*length(G2)),length(G2),length(G2))
    for(i in 11:(10+length(G2))){
      for(j in 11:(10+length(G2))){
        Sigma2[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat2[k,]<-mvrnorm(n = 1, rep(0,length(G2)), Sigma2)
  }
  if (G2.distr=="normal"){
    cov.mat2=cov.mat2
  }else{
    if(G2.distr=="Bernoulli"){
      pvars2 <- pnorm(cov.mat2)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat2 <- qbinom(pvars2, 1, mu.bernoulli)}}
  ###########third group:###################3
  cov.mat3<-matrix(0,n,length(G3))
  for(k in 1:n){
    Sigma3<-matrix(rep(0,length(G3)*length(G3)),length(G3),length(G3))
    for(i in 11:(10+length(G3))){
      for(j in 11:(10+length(G3))){
        Sigma3[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat3[k,]<-mvrnorm(n = 1, rep(0,length(G3)), Sigma3)
  }
  if (G3.distr=="normal"){
    cov.mat3=cov.mat3
  }else{
    if(G3.distr=="Bernoulli"){
      
      pvars3 <- pnorm(cov.mat3)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat3 <- qbinom(pvars3, 1, mu.bernoulli)}}
  #############fourth froup:###############3
  cov.mat4<-matrix(0,n,length(G4))
  for(k in 1:n){
    Sigma4<-matrix(rep(0,length(G4)*length(G4)),length(G4),length(G4))
    for(i in 11:(10+length(G4))){
      for(j in 11:(10+length(G4))){
        Sigma4[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat4[k,]<-mvrnorm(n = 1, rep(0,length(G4)), Sigma4)
  }
  if (G4.distr=="normal"){
    cov.mat4=cov.mat4
  }else{
    if(G4.distr=="Bernoulli"){
      pvars4 <- pnorm(cov.mat4)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat4 <- qbinom(pvars4, 1, mu.bernoulli)}}
  #standardizing Bernoulli ones:
  
  cov.mat=cbind(cov.mat1,cov.mat2,cov.mat3,cov.mat4)
  mycovmean=apply(cov.mat,2,mean)
  p=ncol(cov.mat)
  
  # apply(cov.mat,2,mean)
  #standardize cov.mat:
  cov.mat=scale(cov.mat)
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at age 65
  cens=c(c1,c2)
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:15){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  
genWB.simData.noseed.allcontin.grp.12cov<- function(G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n = n, rho=rho,frailty = T,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  
  
  
  
  cov.mat1<-matrix(0,n,length(G1))
  for(k in 1:n){
    Sigma1<-matrix(rep(0,length(G1)*length(G1)),length(G1),length(G1))
    for(i in 11:(10+length(G1))){
      for(j in 11:(10+length(G1))){
        Sigma1[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat1[k,]<-mvrnorm(n = 1, rep(0,length(G1)), Sigma1)
  }
  if (G1.distr=="normal"){
    cov.mat1=cov.mat1
  }else{
    if(G1.distr=="Bernoulli"){
      pvars1 <- pnorm(cov.mat1)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat1 <- qbinom(pvars1, 1, mu.bernoulli)}}
  #############second group:###########
  cov.mat2<-matrix(0,n,length(G2))
  for(k in 1:n){
    Sigma2<-matrix(rep(0,length(G2)*length(G2)),length(G2),length(G2))
    for(i in 11:(10+length(G2))){
      for(j in 11:(10+length(G2))){
        Sigma2[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat2[k,]<-mvrnorm(n = 1, rep(0,length(G2)), Sigma2)
  }
  if (G2.distr=="normal"){
    cov.mat2=cov.mat2
  }else{
    if(G2.distr=="Bernoulli"){
      pvars2 <- pnorm(cov.mat2)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat2 <- qbinom(pvars2, 1, mu.bernoulli)}}
  ###########third group:###################3
  cov.mat3<-matrix(0,n,length(G3))
  for(k in 1:n){
    Sigma3<-matrix(rep(0,length(G3)*length(G3)),length(G3),length(G3))
    for(i in 11:(10+length(G3))){
      for(j in 11:(10+length(G3))){
        Sigma3[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat3[k,]<-mvrnorm(n = 1, rep(0,length(G3)), Sigma3)
  }
  if (G3.distr=="normal"){
    cov.mat3=cov.mat3
  }else{
    if(G3.distr=="Bernoulli"){
      
      pvars3 <- pnorm(cov.mat3)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat3 <- qbinom(pvars3, 1, mu.bernoulli)}}
  #############fourth froup:###############3
  cov.mat4<-matrix(0,n,length(G4))
  for(k in 1:n){
    Sigma4<-matrix(rep(0,length(G4)*length(G4)),length(G4),length(G4))
    for(i in 11:(10+length(G4))){
      for(j in 11:(10+length(G4))){
        Sigma4[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat4[k,]<-mvrnorm(n = 1, rep(0,length(G4)), Sigma4)
  }
  if (G4.distr=="normal"){
    cov.mat4=cov.mat4
  }else{
    if(G4.distr=="Bernoulli"){
      pvars4 <- pnorm(cov.mat4)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat4 <- qbinom(pvars4, 1, mu.bernoulli)}}
  #standardizing Bernoulli ones:
  
  cov.mat=cbind(cov.mat1,cov.mat2,cov.mat3,cov.mat4)
  mycovmean=apply(cov.mat,2,mean)
  p=ncol(cov.mat)
  
  # apply(cov.mat,2,mean)
  #standardize cov.mat:
  cov.mat=scale(cov.mat)
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  
  data0=data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth)
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at age 65
  cens=c(c1,c2)
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:17){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  
genWB.simData.noseed.allcontin.grp.15cov<- function(G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n = n, rho=rho,frailty = T,risks.status){
  #frailty parameter:
  theta.true = ifelse(frailty == T, 0.25, 0)
  #Weibull parameters:
  # log.kappa1.true = -9.98
  # log.alpha1.true = 1.05
  # log.kappa2.true = -10.01
  # log.alpha2.true = 1.15
  # log.kappa3.true = -5.92
  # log.alpha3.true = 0.92
  # log.kappa1.true =-7.129885398
  # log.alpha1.true =-0.323245915
  # log.kappa2.true = -10.642380920
  # log.alpha2.true = 0.058558940
  # log.kappa3.true = -8.601033114
  # log.alpha3.true = 0.039590459
  log.kappa1.true =weibull.param.log[1]
  log.alpha1.true =weibull.param.log[2]
  log.kappa2.true =weibull.param.log[3]
  log.alpha2.true = weibull.param.log[4]
  log.kappa3.true = weibull.param.log[5]
  log.alpha3.true = weibull.param.log[6]
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  #true beta values:
  # beta1.true = c(-0.8,1,1,0.9,0.8,1,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,-0.8,1,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0.9,0.8,0,0,0,0)
  
  #It is working, but to try with less number of nonzero ones Im gonna try another set of true values below of it.
  beta1.true = beta1.true
  beta2.true = beta2.true
  beta3.true = beta3.true
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  
  
  cens = c(c1,c2)
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  
  
  
  
  cov.mat1<-matrix(0,n,length(G1))
  for(k in 1:n){
    Sigma1<-matrix(rep(0,length(G1)*length(G1)),length(G1),length(G1))
    for(i in 11:(10+length(G1))){
      for(j in 11:(10+length(G1))){
        Sigma1[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat1[k,]<-mvrnorm(n = 1, rep(0,length(G1)), Sigma1)
  }
  if (G1.distr=="normal"){
    cov.mat1=cov.mat1
  }else{
    if(G1.distr=="Bernoulli"){
      pvars1 <- pnorm(cov.mat1)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat1 <- qbinom(pvars1, 1, mu.bernoulli)}}
  #############second group:###########
  cov.mat2<-matrix(0,n,length(G2))
  for(k in 1:n){
    Sigma2<-matrix(rep(0,length(G2)*length(G2)),length(G2),length(G2))
    for(i in 11:(10+length(G2))){
      for(j in 11:(10+length(G2))){
        Sigma2[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat2[k,]<-mvrnorm(n = 1, rep(0,length(G2)), Sigma2)
  }
  if (G2.distr=="normal"){
    cov.mat2=cov.mat2
  }else{
    if(G2.distr=="Bernoulli"){
      pvars2 <- pnorm(cov.mat2)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat2 <- qbinom(pvars2, 1, mu.bernoulli)}}
  ###########third group:###################3
  cov.mat3<-matrix(0,n,length(G3))
  for(k in 1:n){
    Sigma3<-matrix(rep(0,length(G3)*length(G3)),length(G3),length(G3))
    for(i in 11:(10+length(G3))){
      for(j in 11:(10+length(G3))){
        Sigma3[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat3[k,]<-mvrnorm(n = 1, rep(0,length(G3)), Sigma3)
  }
  if (G3.distr=="normal"){
    cov.mat3=cov.mat3
  }else{
    if(G3.distr=="Bernoulli"){
      
      pvars3 <- pnorm(cov.mat3)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat3 <- qbinom(pvars3, 1, mu.bernoulli)}}
  #############fourth froup:###############3
  cov.mat4<-matrix(0,n,length(G4))
  for(k in 1:n){
    Sigma4<-matrix(rep(0,length(G4)*length(G4)),length(G4),length(G4))
    for(i in 11:(10+length(G4))){
      for(j in 11:(10+length(G4))){
        Sigma4[i-10,j-10]<-rho^(abs(i-j))
      }
    }
    cov.mat4[k,]<-mvrnorm(n = 1, rep(0,length(G4)), Sigma4)
  }
  if (G4.distr=="normal"){
    cov.mat4=cov.mat4
  }else{
    if(G4.distr=="Bernoulli"){
      pvars4 <- pnorm(cov.mat4)
      # cov(pvars2);cor(pvars2) We can see that the cov and cor of pvars are very close to our desired Sigma (/rho in Zhaoetal paper)
      cov.mat4 <- qbinom(pvars4, 1, mu.bernoulli)}}
  #standardizing Bernoulli ones:
  
  cov.mat=cbind(cov.mat1,cov.mat2,cov.mat3,cov.mat4)
  mycovmean=apply(cov.mat,2,mean)
  p=ncol(cov.mat)
  
  # apply(cov.mat,2,mean)
  #standardize cov.mat:
  cov.mat=scale(cov.mat)
  # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
  first=as.vector(cov.mat[,1])
  second=as.vector(cov.mat[,2])
  third=as.vector(cov.mat[,3])
  fourth=as.vector(cov.mat[,4])
  fifth=as.vector(cov.mat[,5])
  sixth=as.vector(cov.mat[,6])
  seventh=as.vector(cov.mat[,7])
  eighth=as.vector(cov.mat[,8])
  ninth=as.vector(cov.mat[,9])
  tenth=as.vector(cov.mat[,10])
  eleventh=as.vector(cov.mat[,11])
  twelfth=as.vector(cov.mat[,12])
  thirteenth=as.vector(cov.mat[,13])
  fourteenth=as.vector(cov.mat[,14])
  fifteenth=as.vector(cov.mat[,15])
  
  
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth,fifteenth)
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at age 65
  cens=c(c1,c2)
  Y.tmp = data.frame(simID.LT.noseed(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                                     alpha1.true, alpha2.true, alpha3.true,
                                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:20){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              frailty = frailty,
              WBparamNames = WBparamNames,
              WBpara = WBpara)) 
}  

logLike.weibull.SCR.SM.LT <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE)
{
  ##
  kappa1    <- exp(para[1])
  alpha1 <- exp(para[2])
  kappa2    <- exp(para[3])
  alpha2 <- exp(para[4])
  kappa3    <- exp(para[5])
  alpha3 <- exp(para[6])
  if(frailty == TRUE){
    theta    <- exp(para[7])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, 7, 6)
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
  eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
  eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1 & l < y1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1 & l < y1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0 & l < y1)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0 & l < y1)
  type5 <- as.numeric(delta1 == 1 & delta2 == 1 & y1 <= l & l < y2)
  type6 <- as.numeric(delta1 == 1 & delta2 == 0 & y1 <= l & l < y2)
  ##
  log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  ##
  q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ##
  w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  w.y1.l  <- kappa3*((l-y1)^(alpha3))* exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  k3 <- w.y1.y2 - w.y1.l
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (k1 + k2.y1))))
    logLike2 <- log.h2star.y1 - ((thetaInv + 1) * log(1 + (theta * k2.y1)))  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (k1 + k2.y1))))
    logLike4 <- - thetaInv * log(1 + (theta * k2.y1))  ## Making in terms of y1
    logLike5 <- log.h3star.y2 - ((thetaInv + 1) * log(1 + (theta * k3)))
    logLike6 <- - thetaInv * log(1 + (theta * k3))
  }
  if(frailty == FALSE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
    logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - (k1 + k2.y1)
    logLike4 <- - k2.y1  ## Making in terms of y1
    logLike5 <- log.h3star.y2 - k3
    logLike6 <- - k3
  }
  ##
  loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) + sum(logLike5[type5==1]) + sum(logLike6[type6==1])
  ##
  return(-loglh)
}
logLike.weibull.SCR.SM.LT2 <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE)
{
  ##
  kappa1    <- exp(para[1])
  alpha1 <- exp(para[2])
  kappa2    <- exp(para[3])
  alpha2 <- exp(para[4])
  kappa3    <- exp(para[5])
  alpha3 <- exp(para[6])
  if(frailty == TRUE){
    theta    <- exp(para[7])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, 7, 6)
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
  eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
  eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1 & l < y1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1 & l < y1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0 & l < y1)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0 & l < y1)
  type5 <- as.numeric(delta1 == 1 & delta2 == 1 & y1 <= l & l < y2)
  type6 <- as.numeric(delta1 == 1 & delta2 == 0 & y1 <= l & l < y2)
  ##
  log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  ##
  q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ##
  w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  w.y1.l  <- kappa3*((l-y1)^(alpha3))* exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  k3 <- w.y1.y2 - w.y1.l
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (k1 + k2.y1))))
    logLike2 <- log.h2star.y1 - ((thetaInv + 1) * log(1 + (theta * k2.y1)))  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (k1 + k2.y1))))
    logLike4 <- - thetaInv * log(1 + (theta * k2.y1))  ## Making in terms of y1
    logLike5 <- log.h3star.y2 - ((thetaInv + 1) * log(1 + (theta * k3)))
    logLike6 <- - thetaInv * log(1 + (theta * k3))
  }
  if(frailty == FALSE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
    logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - (k1 + k2.y1)
    logLike4 <- - k2.y1  ## Making in terms of y1
    logLike5 <- log.h3star.y2 - k3
    logLike6 <- - k3
  }
  ##
  
  loglh <- sum((logLike1[type1==1])>-10000) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) + sum(logLike5[type5==1]) + sum(logLike6[type6==1])
  ##
  return(-loglh)
}


# para=para.est
dlogLike.weibull.new <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE)
{
  ##
  kappa1 <- exp(para[1])
  alpha1 <- exp(para[2])
  kappa2 <- exp(para[3])
  alpha2 <- exp(para[4])
  kappa3 <- exp(para[5])
  alpha3 <- exp(para[6])
  if(frailty == TRUE){
    theta    <- exp(para[7])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, 7, 6)
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
  eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
  eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1 & l < y1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1 & l < y1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0 & l < y1)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0 & l < y1)
  type5 <- as.numeric(delta1 == 1 & delta2 == 1 & y1 <= l & l < y2)
  type6 <- as.numeric(delta1 == 1 & delta2 == 0 & y1 <= l & l < y2)
  ##
  log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  ##
  q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ##
  w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  w.y1.l  <- kappa3*(l-y1)^alpha3 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  k3 <- w.y1.y2 - w.y1.l
  ##
  if(frailty == TRUE)
  {
    
    n=NROW(eta.1)
    p=ncol(Xmat1)
    #score_ij-->i(l_i),(j:beta_j)--->score23:dl_2(\phi)/d\beta_3.
    #we have l_i; i=1,2,3,4, and beta_j; j=1,2,3,4.
    A1=(kappa1*(y1)^alpha1*exp(eta.1))-(kappa1*(l)^alpha1*exp(eta.1))
    denom1=1+(theta*(k1+k2.y1))
    denom2=1+(theta*k2.y1)
    A2=(kappa2*(y1)^alpha2 *exp(eta.2))-(kappa2*(l)^alpha2 *exp(eta.2))
    A3=(kappa3*(y2-y1)^alpha3 *exp(eta.3))
    
    
    #each of these below is a n*1 column vector.
    aa1=1-((1+2*theta)*A1/denom1)
    aa2=-(1+theta)*(A1)/(denom2)
    aa3=1-((1+theta)*(A1)/(denom1))
    aa4=-(A1)/(denom2)
    bb1=-(1+2*theta)*(A2)/(denom1)
    bb2=1-(((1+theta)*A2)/(denom2))
    bb3=-(1+theta)*(A2)/(denom1)
    bb4=-(A2)/(denom2)
    cc1=1-(((1+2*theta)*A3)/(denom1))
    cc3=-(1+theta)*(A3)/(denom1)
    
    
    #now, we multiply a X_ik; k=1,2,3, i=1,2,...,n in those aa1, aa2, ... to construct the components of the dloglikelihood function.
    #X_ik is a p-vector in the column form.
    #So, by multiplying, we have a matrix of n*p for each of score_k; k=1,2,3.
    #each score_k; k=1,2,3 consists of 4 parts that are the components in the likelihood function.
    #4 components of score1:
    score11=matrix(NA,n,p)
    score12=matrix(NA,n,p)
    score13=matrix(NA,n,p)
    score14=matrix(NA,n,p)
    #4 components of score2:
    score21=matrix(NA,n,p)
    score22=matrix(NA,n,p)
    score23=matrix(NA,n,p)
    score24=matrix(NA,n,p)
    #2 components of score3 (2 of them are zero):
    score31=matrix(NA,n,p)
    score33=matrix(NA,n,p)
    
    for (i in 1:n){
      #each row in the n*p matrix is the multipication of each element ofaa1, bb1, or etc for i=1,2,...,n in the X_ik
      score11[i,]=as.vector(aa1[i]*Xmat1[i,])
      score12[i,]=as.vector(aa2[i]*Xmat1[i,])
      score13[i,]=as.vector(aa3[i]*Xmat1[i,])
      score14[i,]=as.vector(aa4[i]*Xmat1[i,])
      
      score21[i,]=as.vector(bb1[i]*Xmat2[i,])
      score22[i,]=as.vector(bb2[i]*Xmat2[i,])
      score23[i,]=as.vector(bb3[i]*Xmat2[i,])
      score24[i,]=as.vector(bb4[i]*Xmat2[i,])
      
      score31[i,]=as.vector(cc1[i]*Xmat3[i,])
      score33[i,]=as.vector(cc3[i]*Xmat3[i,])
    }
    #each of the score1, score2, and score3 below are of .
    #Now, by score_k; k=1,2,3, we have 3 matrices each of n*p dimension:
    score1=(score11*(type1==1))+(score12*(type2==1))+(score13*(type3==1))+(score14*(type4==1))
    score2=(score21*(type1==1))+(score22*(type2==1))+(score23*(type3==1))+(score24*(type4==1))
    score3=(score31*(type1==1))+(score33*(type3==1))
    
    #now we sum over i=1,2,...,n which means that we sum over each column so that we have 1*p vector for each scre1, score2, and score3:
    score1.final=colSums(score1)
    score2.final=colSums(score2)
    score3.final=colSums(score3)
    #finally, score is a vecrtor consisiting of 3 p-vectors:
    score=c(score1.final,score2.final,score3.final)
    
  }
  
  
  return(-score)
}


ddlogLike.weibull.new <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE)
{
  ##
  kappa1 <- exp(para[1])
  alpha1 <- exp(para[2])
  kappa2 <- exp(para[3])
  alpha2 <- exp(para[4])
  kappa3 <- exp(para[5])
  alpha3 <- exp(para[6])
  if(frailty == TRUE){
    theta    <- exp(para[7])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, 7, 6)
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
  eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
  eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1 & l < y1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1 & l < y1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0 & l < y1)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0 & l < y1)
  type5 <- as.numeric(delta1 == 1 & delta2 == 1 & y1 <= l & l < y2)
  type6 <- as.numeric(delta1 == 1 & delta2 == 0 & y1 <= l & l < y2)
  ##
  log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  ##
  q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ##
  w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  w.y1.l  <- kappa3*(l-y1)^alpha3 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  k3 <- w.y1.y2 - w.y1.l
  
  
  if(frailty == TRUE)
  {
    n=NROW(eta.1)
    p=ncol(Xmat1)
    #score_ij-->i(l_i),(j:beta_j)--->score23:dl_2(\phi)/d\beta_3.
    #we have l_i; i=1,2,3,4, and beta_j; j=1,2,3,4.
    A1=(theta*kappa1*(y1)^alpha1*exp(eta.1))-(theta*kappa1*(l)^alpha1*exp(eta.1))
    A2=(theta*kappa2*(y1)^alpha2 *exp(eta.2))-(theta*kappa2*(l)^alpha2 *exp(eta.2))
    A3=theta*(kappa3*(y2-y1)^alpha3 *exp(eta.3))
    A4=kappa1*(y1^alpha1-l^alpha1)*exp(eta.1)
    A5=(kappa1*(y1)^alpha1*exp(eta.1))-(kappa1*(l)^alpha1*exp(eta.1))
    A6=(kappa2*(y1)^alpha2 *exp(eta.2))-(kappa2*(l)^alpha2 *exp(eta.2))
    A7=(kappa3*(y2-y1)^alpha3 *exp(eta.3))
    
    denom1=1+(theta*(k1+k2.y1))
    denom2=1+(theta*k2.y1)
    
    #All n*1 column vectors:
    #dbeta1beta1
    B111=-(1+2*theta)*A4*((denom1-A1)/(denom1^2))
    B211=-(theta+1)*A4*((denom2-A1)/(denom2^2))
    B311=-(1+theta)*A4*((denom1-A1)/(denom1^2))
    B411=-A4*((denom2-A1)/(denom2^2))
    #dbeta1beta2:
    B112=(2*theta+1)*((A2*A5)/(denom1^2))
    B212=(theta+1)*((A2*A5)/(denom2^2))
    B312=(theta+1)*((A2*A5)/(denom1^2))
    B412=(A2*A5)/(denom2^2)
    #dbeta1beta3:
    B113=(2*theta+1)*((A3*A4)/(denom1^2))
    B213=0
    B313=(1+theta)*((A3*A4)/(denom1^2))
    B413=0  
    #dbeta2beta3:
    B123=(1+2*theta)*((A3*A6)/(denom1^2))
    B223=0
    B323=(1+theta)*((A3*A6)/(denom1^2))
    B423=0
    #dbeta2beta2
    B122=-(1+2*theta)*((A6*(denom1-A2))/denom1^2)
    B222=-(1+theta)*((A6*(denom2-A2))/denom2^2)
    B322=-(1+theta)*((A6*(denom1-A2))/denom1^2)
    B422=-((A6*(denom2-A2))/denom2^2)
    #dbeta3beta3:
    B133=-(1+2*theta)*((A7*(denom1-A3))/(denom1^2))
    B333=-(1+theta)*((A7*(denom1-A3))/(denom1^2))
    
    
    dscore111=matrix(0,p,p)
    dscore211=matrix(0,p,p)
    dscore311=matrix(0,p,p)
    dscore411=matrix(0,p,p)
    dscore112=matrix(0,p,p)
    dscore212=matrix(0,p,p)
    dscore312=matrix(0,p,p)
    dscore412=matrix(0,p,p)
    dscore113=matrix(0,p,p)
    dscore313=matrix(0,p,p)
    dscore123=matrix(0,p,p)
    dscore323=matrix(0,p,p)
    dscore122=matrix(0,p,p)
    dscore222=matrix(0,p,p)
    dscore322=matrix(0,p,p)
    dscore422=matrix(0,p,p)
    dscore133=matrix(0,p,p)
    dscore333=matrix(0,p,p)
    
    for (i in 1:n){
      #all p*p matrices:
      #dbeta1dbeta1:
      dscore111=dscore111+(B111[i]*(Xmat1[i,]%*%t(Xmat1[i,])))*((type1==1)[i])
      dscore211=dscore211+(B211[i]*(Xmat1[i,]%*%t(Xmat1[i,])))*((type2==1)[i])
      dscore311=dscore311+(B311[i]*(Xmat1[i,]%*%t(Xmat1[i,])))*((type3==1)[i])
      dscore411=dscore411+(B411[i]*(Xmat1[i,]%*%t(Xmat1[i,])))*((type4==1)[i])
      #dbeta1dbeta2:
      dscore112=dscore112+(B112[i]*(Xmat1[i,]%*%t(Xmat2[i,])))*((type1==1)[i])
      dscore212=dscore212+(B212[i]*(Xmat1[i,]%*%t(Xmat2[i,])))*((type2==1)[i])
      dscore312=dscore312+(B312[i]*(Xmat1[i,]%*%t(Xmat2[i,])))*((type3==1)[i])
      dscore412=dscore412+(B412[i]*(Xmat1[i,]%*%t(Xmat2[i,])))*((type4==1)[i])
      #dbeta1dbeta3:
      dscore113=dscore113+(B113[i]*(Xmat1[i,]%*%t(Xmat3[i,])))*((type1==1)[i])
      # dscore213=dscore213+(B213[i]*(Xmat1[i,]%*%t(Xmat3[i,])))*((type2==1)[i])
      dscore313=dscore313+(B313[i]*(Xmat1[i,]%*%t(Xmat3[i,])))*((type3==1)[i])
      # dscore413=dscore413+(B413[i]*(Xmat1[i,]%*%t(Xmat3[i,])))*((type4==1)[i])
      #dbeta2dbeta3:
      dscore123=dscore123+(B123[i]*(Xmat2[i,]%*%t(Xmat3[i,])))*((type1==1)[i])
      # dscore223=dscore223+(B223[i]*(Xmat2[i,]%*%t(Xmat3[i,])))*((type2==1)[i])
      dscore323=dscore323+(B323[i]*(Xmat2[i,]%*%t(Xmat3[i,])))*((type3==1)[i])
      # dscore423=dscore423+(B423[i]*(Xmat2[i,]%*%t(Xmat3[i,])))*((type4==1)[i])
      #dbeta2dbeta2:
      dscore122=dscore122+(B122[i]*(Xmat2[i,]%*%t(Xmat2[i,])))*((type1==1)[i])
      dscore222=dscore222+(B222[i]*(Xmat2[i,]%*%t(Xmat2[i,])))*((type2==1)[i])
      dscore322=dscore322+(B322[i]*(Xmat2[i,]%*%t(Xmat2[i,])))*((type3==1)[i])
      dscore422=dscore422+(B422[i]*(Xmat2[i,]%*%t(Xmat2[i,])))*((type4==1)[i])
      #dbeta3dbeta3:
      dscore133=dscore133+(B133[i]*(Xmat3[i,]%*%t(Xmat3[i,])))*((type1==1)[i])
      # dscore233=dscore233+(B233[i]*(Xmat3[i,]%*%t(Xmat3[i,])))*((type2==1)[i])
      dscore333=dscore333+(B333[i]*(Xmat3[i,]%*%t(Xmat3[i,])))*((type3==1)[i])
      # dscore433=dscore433+(B433[i]*(Xmat3[i,]%*%t(Xmat3[i,])))*((type4==1)[i])
    }
    dscore11=dscore111+dscore211+dscore311+dscore411
    dscore12=dscore112+dscore212+dscore312+dscore412
    dscore13=dscore113+dscore313
    dscore21=t(dscore12)
    dscore22=dscore122+dscore222+dscore322+dscore422
    dscore23=dscore123+dscore323
    dscore31=t(dscore13)
    dscore32=t(dscore23)
    dscore33=dscore133+dscore333
  }
  
  ##
  #The score function or the first derivative of the loglikelihood function:
  dscore=rbind(cbind(dscore11,dscore12,dscore13),cbind(dscore21,dscore22,dscore23),cbind(dscore31,dscore32,dscore33))
  ##
  return(-dscore)
}
## Fit model: illness-death, shared frailty, left-truncated data
## Weibull baseline hazards
FreqID.LT <- function(Y, lin.pred, data, model = "semi-Markov", startVals, frailty=TRUE, method)
{	
  
  ##
  y1     <- as.vector(Y[,1])
  delta1 <- as.vector(Y[,2])
  y2     <- as.vector(Y[,3])
  delta2 <- as.vector(Y[,4])
  l      <- as.vector(Y[,5])
  # Xmat1=as.matrix(data[,(6:(p+5))])
  # Xmat2=as.matrix(data[,(6:(p+5))])
  # Xmat3=as.matrix(data[,(6:(p+5))])
  Xmat1  <- as.matrix(model.frame(lin.pred[[1]], data=data))
  Xmat2  <- as.matrix(model.frame(lin.pred[[2]], data=data))
  Xmat3  <- as.matrix(model.frame(lin.pred[[3]], data=data))
  ##
  
  fit.survreg.1 <- survreg(as.formula(paste("Surv(y1, delta1) ", as.character(lin.pred[[1]])[1], as.character(lin.pred[[1]])[2])), dist="weibull", data=data,control = list(maxiter=90))
  
  fit.survreg.2 <- survreg(as.formula(paste("Surv(y2, delta2) ", as.character(lin.pred[[2]])[1], as.character(lin.pred[[2]])[2])), dist="weibull", data=data,control = list(maxiter=90))
  data.delta1_1 = data[delta1==1,]
  data.delta1_1$y2.m.y1 = y2[delta1==1] - y1[delta1==1]
  fit.survreg.3 <- survreg(as.formula(paste("Surv(y2.m.y1, delta2) ", as.character(lin.pred[[3]])[1], as.character(lin.pred[[3]])[2])), dist="weibull", data=data.delta1_1,control = list(maxiter=90))
  alpha1      <- 1 / fit.survreg.1$scale
  alpha2      <- 1 / fit.survreg.2$scale
  alpha3     	<- 1 / fit.survreg.3$scale
  
  print(alpha1)
  print(alpha2)
  print(alpha3)
  
  print(coef(fit.survreg.1)[1])
  print(coef(fit.survreg.2)[1])
  print(coef(fit.survreg.3)[1])
  if (is.null(startVals)==T){
    startVals     <- c(-alpha1*coef(fit.survreg.1)[1], log(alpha1),
                       -alpha2*coef(fit.survreg.2)[1], log(alpha2),
                       -alpha3*coef(fit.survreg.3)[1], log(alpha3))
    if(frailty == TRUE) startVals <- c(startVals, 0.5)
    startVals     <- c(startVals,
                       -coef(fit.survreg.1)[-1] * alpha1,
                       -coef(fit.survreg.2)[-1] * alpha2,
                       -coef(fit.survreg.3)[-1] * alpha3)
  }
  
  if(model == "semi-Markov")
  {
    if (method == "optim"){
      cat("Fitting illness-death model with Weibull baseline hazards ... this should take < 1 min \n")
      logLike <- function(p) logLike.weibull.SCR.SM.LT(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                       Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty)
      optim.control = list(REPORT = 50)
      fit1 <- optim(startVals, #* runif(length(startVals), 0.9, 1.1), 
                    logLike, hessian = TRUE, method="Nelder-Mead", control = optim.control)
      value <- list(estimate=fit1$par, H=fit1$hessian, logLike=-fit1$value, code=fit1$convergence)#, Xmat=list(Xmat1, Xmat2, Xmat3))
    }
    if (method == "nlm"){
      cat("Fitting illness-death model with Weibull baseline hazards ... this should take < 1 min \n")
      fit1  <- suppressWarnings(nlm(logLike.weibull.SCR.SM.LT, p=startVals,
                                    y1=y1, delta1=delta1, y2=y2, delta2=delta2, Xmat1=as.matrix(Xmat1), Xmat2=as.matrix(Xmat2), Xmat3=as.matrix(Xmat3),
                                    l = l, frailty=frailty,
                                    iterlim=1000, hessian=TRUE))
      value <- list(estimate=fit1$est, H=fit1$hessian, logLike=-fit1$minimum, code=fit1$code)
    }
  }
  ##
  if(model == "semi-Markov")
  {
    class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
  }
  
  return(value)
  ##
  invisible()
}
FreqID.LT.oracle.overlapping <- function(Y, lin.pred.oracle, data.oracle.1,data.oracle.2,data.oracle.3, model = "semi-Markov", startVals, frailty, method)
{	
  
  ##
  y1     <- as.vector(Y[,1])
  delta1 <- as.vector(Y[,2])
  y2     <- as.vector(Y[,3])
  delta2 <- as.vector(Y[,4])
  l      <- as.vector(Y[,5])
  # Xmat1=as.matrix(data[,(6:(p+5))])
  # Xmat2=as.matrix(data[,(6:(p+5))])
  # Xmat3=as.matrix(data[,(6:(p+5))])
  
  
  Xmat1  <- as.matrix(model.frame(lin.pred.oracle[[1]], data=data.oracle.1))
  Xmat2  <- as.matrix(model.frame(lin.pred.oracle[[2]], data=data.oracle.2))
  Xmat3  <- as.matrix(model.frame(lin.pred.oracle[[3]], data=data.oracle.3))
  
  
  
  ##
  
  
  fit.survreg.1<- survreg(as.formula(paste("Surv(y1, delta1) ", as.character(lin.pred.oracle[[1]])[1], as.character(lin.pred.oracle[[1]])[2])), dist="weibull", data=data.oracle.1,control = list(maxiter=90))
  fit.survreg.2<- survreg(as.formula(paste("Surv(y2, delta2) ", as.character(lin.pred.oracle[[2]])[1], as.character(lin.pred.oracle[[2]])[2])), dist="weibull", data=data.oracle.2,control = list(maxiter=90))
  data.delta1_1 = as.data.frame(as.matrix(cbind(((data.oracle.1[delta1==1,])[,(1:5)]),(Xmat3[delta1==1,]))))
  data.delta1_1$y2.m.y1 = y2[delta1==1] - y1[delta1==1]
  fit.survreg.3 <- survreg(as.formula(paste("Surv(y2.m.y1, delta2) ", as.character(lin.pred.oracle[[3]])[1], as.character(lin.pred.oracle[[3]])[2])), dist="weibull", data=data.delta1_1,control = list(maxiter=90) )
  alpha1      <- 1 / fit.survreg.1$scale
  alpha2      <- 1 / fit.survreg.2$scale
  alpha3     	<- 1 / fit.survreg.3$scale
  
  # 
  # print(alpha1)
  # print(alpha2)
  # print(alpha3)
  # 
  # print(coef(fit.survreg.1)[1])
  # print(coef(fit.survreg.2)[1])
  # print(coef(fit.survreg.3)[1])
  if (is.null(startVals)==T){
    startVals     <- c(-alpha1*coef(fit.survreg.1)[1], log(alpha1),
                       -alpha2*coef(fit.survreg.2)[1], log(alpha2),
                       -alpha3*coef(fit.survreg.3)[1], log(alpha3))
    if(frailty == TRUE) startVals <- c(startVals, 0.5)
    startVals     <- c(startVals,
                       -coef(fit.survreg.1)[-1] * alpha1,
                       -coef(fit.survreg.2)[-1] * alpha2,
                       -coef(fit.survreg.3)[-1] * alpha3)
  }
  
  if(model == "semi-Markov")
  {
    if (method == "optim"){
      cat("Fitting illness-death model with Weibull baseline hazards ... this should take < 1 min \n")
      logLike <- function(p) logLike.weibull.SCR.SM.LT(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                       Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty)
      optim.control = list(REPORT = 50)
      fit1 <- optim(startVals, #* runif(length(startVals), 0.9, 1.1), 
                    logLike, hessian = TRUE, method="Nelder-Mead", control = optim.control)
      value <- list(estimate=fit1$par, H=fit1$hessian, logLike=-fit1$value, code=fit1$convergence)#, Xmat=list(Xmat1, Xmat2, Xmat3))
    }
    if (method == "nlm"){
      cat("Fitting illness-death model with Weibull baseline hazards ... this should take < 1 min \n")
      fit1  <- suppressWarnings(nlm(logLike.weibull.SCR.SM.LT, p=startVals,
                                    y1=y1, delta1=delta1, y2=y2, delta2=delta2, Xmat1=as.matrix(Xmat1), Xmat2=as.matrix(Xmat2), Xmat3=as.matrix(Xmat3),
                                    l = l, frailty=frailty,
                                    iterlim=1000, hessian=TRUE))
      value <- list(estimate=fit1$est, H=fit1$hessian, logLike=-fit1$minimum, code=fit1$code)
    }
  }
  ##
  if(model == "semi-Markov")
  {
    class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
  }
  
  return(value)
  ##
  invisible()
}
