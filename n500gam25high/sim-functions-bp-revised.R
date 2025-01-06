library(dplyr)
library(survival)
require(stats)
library(splines2)
library(pracma)  ## for numerical differentiation
library(MASS)
# Data Generation ---------------------------------------------------------


## For estimating B-spline basis coefficients
obj.fn.1 <- function(eta.vec, spline.pred, h0.truth){ 
  ##
  spline.vals <- sapply(1:nrow(spline.pred), function(i) sum(spline.pred[i,]*(eta.vec)))
  sq.diffs <- (log(h0.truth)-spline.vals)^2
  return(sum(sq.diffs))
}

## Simulate data from illness-death, shared frailty
## Adapted from simID in SemiCompRisks package
## lt.type: draw left-truncation times from "unif" (uniform) or "norm" (normal)
##          - if "unif", lt=c(a, b) so that L ~ unif(n, a, b)
##          - if "norm", lt=c(mean, sd), so that L ~ Norm(mean, sd)
simID.LT <- function(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                     alpha1.true, alpha2.true, alpha3.true,
                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt.type = "unif", lt,risks.status,frailty_distribution)
{
  
  n <- dim(x1)[1]
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]
  p3 <- dim(x3)[2]
  
  
  if (frailty_distribution=="LogNormal"){
    gamma.true=rlnorm(n,0,theta.true)
  }
  if(frailty_distribution=="Gamma"){
    if(theta.true >0)
    {
      gamma.true <- rgamma(n, 1/theta.true, 1/theta.true)
    }
    if(theta.true == 0)
    {
      gamma.true <- rep(1, n)
    }
  }
  # if(theta.true >0)
  # {
  #   gamma.true <- rgamma(n, 1/theta.true, 1/theta.true)
  # }
  # if(theta.true == 0)
  # {
  #   gamma.true <- rep(1, n)
  # }
  
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

## For plotting WB baseline hazard function
WB.haz <- function(log.kappa, log.alpha, t){
  kappa = exp(log.kappa); alpha = exp(log.alpha)
  h = alpha * kappa * t^(alpha - 1)
  return(h)
}

## Generating simulated data used in simulation study

# genWB.simData.10cov.BP.BS <- function(weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution){
#   
#   # theta.true = ifelse(frailty == T, 0.25, 0)
#   theta.true = ifelse(frailty == T, theta.given, 0)
#   
#   log.kappa1.true =weibull.param.log[1]
#   log.alpha1.true =weibull.param.log[2]
#   log.kappa2.true =weibull.param.log[3]
#   log.alpha2.true = weibull.param.log[4]
#   log.kappa3.true = weibull.param.log[5]
#   log.alpha3.true = weibull.param.log[6]  
#   ##
#   kappa1.true = exp(log.kappa1.true)
#   alpha1.true = exp(log.alpha1.true)
#   kappa2.true = exp(log.kappa2.true)
#   alpha2.true = exp(log.alpha2.true)
#   kappa3.true = exp(log.kappa3.true)
#   alpha3.true = exp(log.alpha3.true)
#   
#   # beta1.true = c(-0.03021803)  
#   # beta2.true = c(-0.33064510)  
#   # beta3.true = c(-0.10652843)
#   
#   beta1.true=beta1.true
#   beta2.true=beta2.true
#   beta3.true=beta3.true
#   
#   
#   cens = c(c1,c2)
#   
#   
#   
#   m1=m[1]
#   m2=m[2]
#   m3=m[3]
#   
#   ## Male 
#   
#   
#   ## B-spline baseline hazard specifications
#   n.internalKnots.1 <- 1
#   Bspline.degree.1 <- 1
#   num.Bspline.params.1 <- (n.internalKnots.1) + (Bspline.degree.1 + 1)
#   
#   n.internalKnots.2 <- 1
#   Bspline.degree.2 <- 1
#   num.Bspline.params.2 <- (n.internalKnots.2) + (Bspline.degree.2 + 1)
#   
#   n.internalKnots.3 <- 1  ## Trying to add more flexibility - the previously estimated BH was off!
#   Bspline.degree.3 <- 2
#   num.Bspline.params.3 <- (n.internalKnots.3) + (Bspline.degree.3 + 1)
#   
#   lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth))
#   if (risks.status=="fullyshared"){
#     lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
#   }
#   if (risks.status=="overlapping"){
#     lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
#   }
#   
#   num.Bspline.params.1 = (n.internalKnots.1) + (Bspline.degree.1 + 1)
#   num.Bspline.params.2 = (n.internalKnots.2) + (Bspline.degree.2 + 1)
#   num.Bspline.params.3 = (n.internalKnots.3) + (Bspline.degree.3 + 1)
#   
#   n.params.BS <- ifelse(frailty == T, num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true) + 1, 
#                         num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true))
#   
#   BSparamNames <- c(paste0("eta1.", c(0:(m1))), paste0("eta2.", c(0:(m2))), paste0("eta3.", c(0:(m3))))
#   BPparamNames=c(paste0("phi1.", c(0:(num.Bspline.params.1-1))), paste0("phi2.", c(0:(num.Bspline.params.2-1))), paste0("phi3.", c(0:(num.Bspline.params.3-1))))
#   WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
#   if (frailty == T){
#     BSparamNames <- c(BSparamNames, "log(theta)")
#     WBparamNames <- c(WBparamNames, "log(theta)") 
#     BPparamNames=c(BPparamNames, "log(theta)")
#   }
#   BSparamNames <- c(BSparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
#   BPparamNames <- c(BPparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
#   WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
#   
#   ## DATA
#   
#   library(MASS)
#   cov.mat<-matrix(0,n,10)
#   for(k in 1:n){
#     Sigma<-matrix(rep(0,10*10),10,10)
#     for(i in 11:20){
#       for(j in 11:20){
#         Sigma[i-10,j-10]<-rho^(abs(i-j))
#       }
#     }
#     cov.mat[k,]<-mvrnorm(n = 1, rep(0, 10), Sigma)
#   }
#   # cov.mat=(cov.mat-min(cov.mat))/(max(cov.mat)-min(cov.mat))
#   first=as.vector(cov.mat[,1])
#   second=as.vector(cov.mat[,2])
#   third=as.vector(cov.mat[,3])
#   fourth=as.vector(cov.mat[,4])
#   fifth=as.vector(cov.mat[,5])
#   sixth=as.vector(cov.mat[,6])
#   seventh=as.vector(cov.mat[,7])
#   eighth=as.vector(cov.mat[,8])
#   ninth=as.vector(cov.mat[,9])
#   tenth=as.vector(cov.mat[,10])
#   
#   data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth)
#   
#   x1 = x2 = x3 = as.matrix(data0)
#   
#   lt = c(lt1,lt2) # Origin starts at study entry
#   
#   Y.tmp = data.frame(simID.LT(x1, x2, x3, beta1.true, beta2.true, beta3.true,
#                               alpha1.true, alpha2.true, alpha3.true,
#                               kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt,frailty_distribution), data0)
#   ## Make sure left-truncated data
#   Y = Y.tmp %>% filter(y1 > L)
#   for (i in 6:15){
#     Y[,i]=scale(Y[,i])
#   }
#   
#   data = data.frame(Y)
#   
#   delta1 = Y$delta1; delta2 = Y$delta2
#   y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
#   n=length(y1)
#   ## Estimate B-spline params
#   ## log h0(t) = phi0 * B0(t) + ... + phik * Bk(t)
#   bdy.knots.b.1 = c(0, max(y1)) 
#   bdy.knots.b.2 = c(0, max(y2))
#   bdy.knots.b.3.y2my1 = c(0, max(y2-y1))
#   
#   # #1 tarnsition:
#   h0.1.truth.event <- alpha1.true*kappa1.true*y1[delta1==1]^(alpha1.true-1)
#   
#   ## 1-transition-BS:
#   knot.loc.1 <- quantile(y1[delta1==1], ppoints(n.internalKnots.1))
#   b.1.event.bs <- bSpline(y1[delta1==1], knots = knot.loc.1, degree = Bspline.degree.1, intercept=TRUE, Boundary.knots=bdy.knots.b.1)
#   b.1.bs <- predict(b.1.event.bs, y1)
#   phi.1.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.1), spline.pred=b.1.event.bs, h0.truth=h0.1.truth.event)$estimate
#   
#   #1-transition-BP:
#   b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
#   b.1.bp=predict(b.1.event.bp,y1)
#   # b.1.bp=BP.basis.finder(n,m1,y1)
#   # n.event.1=length(y1[delta1==1])
#   # b.1.bp.event=BP.basis.finder(n.event.1,m1,y1[delta1==1])
#   phi.1.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m1+1)),spline.pred=b.1.event.bp,h0.truth=h0.1.truth.event)$estimate
#   start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength, method="nls")
#   
#   # #2 transition:
#   h0.2.truth.event <- alpha2.true*kappa2.true*y2[delta2==1]^(alpha2.true-1)
#   ## 2-transition-BS:
#   knot.loc.2 <- quantile(y2[delta2==1], ppoints(n.internalKnots.2))
#   b.2.event.bs <- bSpline(y2[delta2==1], knots = knot.loc.2, degree = Bspline.degree.2, intercept=TRUE, Boundary.knots=bdy.knots.b.2)
#   b.2.bs <- predict(b.2.event.bs, y2)
#   phi.2.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.2), spline.pred=b.2.event.bs, h0.truth=h0.2.truth.event)$estimate
#   #2-transition-BP:
#   b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
#   b.2.bp=predict(
#     b.2.event.bp,y2)
#   # b.2.bp=BP.basis.finder(n,m2,y2)
#   # n.event.2=length(y2[delta2==1])
#   # b.2.bp.event=BP.basis.finder(n.event.2,m1,y2[delta2==1])
#   phi.2.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m2+1)),spline.pred=b.2.event.bp,h0.truth=h0.2.truth.event)$estimate
#   start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength, method="nls")
#   
#   ## 3-transition:
#   h0.3.truth.event <- alpha3.true*kappa3.true*(y2-y1)[delta1==1]^(alpha3.true-1)
#   
#   # #3 transition-BS:
#   knot.loc.3 <- quantile((y2-y1)[delta1==1], ppoints(n.internalKnots.3))
#   b.3.event.bs <- bSpline((y2-y1)[delta1==1], knots = knot.loc.3, degree = Bspline.degree.3, intercept=TRUE, Boundary.knots=bdy.knots.b.3.y2my1)
#   b.3.y2my1.bs <- predict(b.3.event.bs, y2-y1)
#   phi.3.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.3), spline.pred=b.3.event.bs, h0.truth=h0.3.truth.event)$estimate
#   
#   #3-transition-BP:
#   # b.3.bp=BP.basis.finder(n,m3,(y2-y1))
#   # n.event.3=length((y2-y1)[delta1==1])
#   # b.3.bp.event=BP.basis.finder(n.event.3,m3,(y2-y1)[delta1==1])
#   b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
#   b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
#   phi.3.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m3+1)),spline.pred=b.3.event.bp,h0.truth=h0.3.truth.event)$estimate
#   start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength, method="nls")
#   
#   ## Start values for b-spline simulations    
#   startVals.bs = c(phi.1.truth.bs,
#                    phi.2.truth.bs,
#                    phi.3.truth.bs)
#   if (frailty == T) startVals.bs = c(startVals.bs, log(theta.true))
#   startVals.bs = c(startVals.bs, beta1.true, beta2.true, beta3.true)
#   
#   
#   beta1.start.bp=start1.bp[(length(phi.1.truth.bp)+1):length(start1.bp)]
#   beta2.start.bp=start2.bp[(length(phi.2.truth.bp)+1):length(start2.bp)]
#   beta3.start.bp=start3.bp[(length(phi.3.truth.bp)+1):length(start3.bp)]
#   
#   ## Start values for Bernstein simulations    
#   startVals.bp = c(phi.1.truth.bp,
#                    phi.2.truth.bp,
#                    phi.3.truth.bp)
#   if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true)+rnorm(1, mean=0, sd=0.1))
#   # if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true))
#   startVals.bp = c(startVals.bp, beta1.start.bp, beta2.start.bp, beta3.start.bp)
#   
#   
#   phi.truth.BS=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
#   phi.truth.BP=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
#   
#   #Weibull:
#   WBpara = c(log(kappa1.true), log(alpha1.true),
#              log(kappa2.true), log(alpha2.true),
#              log(kappa3.true), log(alpha3.true)) 
#   if (frailty == T) WBpara = c(WBpara, log(theta.true))
#   WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
#   
#   #BS:
#   BSpara=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
#   if (frailty == T) BSpara=c(BSpara,log(theta.true))
#   BSpara=c(BSpara,beta1.true,beta2.true,beta3.true)
#   
#   #BP:
#   BPpara=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
#   if (frailty == T) BPpara=c(BPpara,log(theta.true))
#   BPpara=c(BPpara,beta1.true,beta2.true,beta3.true)
#   
#   
#   return(list(Y = Y,
#               lin.pred = lin.pred,
#               lin.pred.oracle=lin.pred.oracle,
#               data = data,
#               startVals.bs = startVals.bs,
#               startVals.bp=startVals.bp,
#               frailty = frailty,
#               b.1.bs=b.1.bs, b.2.bs=b.2.bs, b.3.y2my1.bs=b.3.y2my1.bs,
#               bdy.knots.b.1=bdy.knots.b.1, 
#               bdy.knots.b.2=bdy.knots.b.2, 
#               bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1, 
#               WBparamNames = WBparamNames,
#               BSparamNames = BSparamNames,
#               BPparamNames = BPparamNames,
#               num.Bspline.params.1 = num.Bspline.params.1,
#               num.Bspline.params.2 = num.Bspline.params.2,
#               num.Bspline.params.3 = num.Bspline.params.3,
#               m=m,
#               WBpara = WBpara, 
#               BSpara=BSpara,
#               BPpara=BPpara,
#               b.1.bp=b.1.bp,b.2.bp=b.2.bp,b.3.y2my1.bp=b.3.y2my1.bp,
#               phi.truth.BS=phi.truth.BS,phi.truth.BP=phi.truth.BP))
# }  
genWB.simData.12cov.BP.BS<- function(weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution){
  #frailty parameter:
  # theta.true = ifelse(frailty == T, 0.25, 0)
  theta.true = ifelse(frailty == T, theta.given, 0)
  #Weibull parameters:
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
  
  
  
  m1=m[1]
  m2=m[2]
  m3=m[3]
  
  n.internalKnots.1 <- 1
  Bspline.degree.1 <- 1
  num.Bspline.params.1 <- (n.internalKnots.1) + (Bspline.degree.1 + 1)
  
  n.internalKnots.2 <- 1
  Bspline.degree.2 <- 1
  num.Bspline.params.2 <- (n.internalKnots.2) + (Bspline.degree.2 + 1)
  
  n.internalKnots.3 <- 1  ## Trying to add more flexibility - the previously estimated BH was off!
  Bspline.degree.3 <- 2
  num.Bspline.params.3 <- (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth))
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  
  
  num.Bspline.params.1 = (n.internalKnots.1) + (Bspline.degree.1 + 1)
  num.Bspline.params.2 = (n.internalKnots.2) + (Bspline.degree.2 + 1)
  num.Bspline.params.3 = (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  n.params.BS <- ifelse(frailty == T, num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true) + 1, 
                        num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true))
  
  BSparamNames <- c(paste0("eta1.", c(0:(m1))), paste0("eta2.", c(0:(m2))), paste0("eta3.", c(0:(m3))))
  BPparamNames=c(paste0("phi1.", c(0:(num.Bspline.params.1-1))), paste0("phi2.", c(0:(num.Bspline.params.2-1))), paste0("phi3.", c(0:(num.Bspline.params.3-1))))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    BSparamNames <- c(BSparamNames, "log(theta)")
    WBparamNames <- c(WBparamNames, "log(theta)") 
    BPparamNames=c(BPparamNames, "log(theta)")
  }
  BSparamNames <- c(BSparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  BPparamNames <- c(BPparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
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
    cov.mat[k,]<-mvrnorm(n = 1, rep(0,12), Sigma)
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
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  
  Y.tmp = data.frame(simID.LT(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                              alpha1.true, alpha2.true, alpha3.true,
                              kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens,lt.type="unif", lt = lt,risks.status,frailty_distribution), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:17){
    Y[,i]=scale(Y[,i])
  }
  
  # Y.oracle=Y[,(1:10)]
  # data.oracle=data.frame(Y.oracle)
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  n=length(y1)
  bdy.knots.b.1 = c(0, max(y1)) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2-y1))
  
  # #1 tarnsition:
  h0.1.truth.event <- alpha1.true*kappa1.true*y1[delta1==1]^(alpha1.true-1)
  ## 1-transition-BS:
  knot.loc.1 <- quantile(y1[delta1==1], ppoints(n.internalKnots.1))
  b.1.event.bs <- bSpline(y1[delta1==1], knots = knot.loc.1, degree = Bspline.degree.1, intercept=TRUE, Boundary.knots=bdy.knots.b.1)
  b.1.bs <- predict(b.1.event.bs, y1)
  phi.1.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.1), spline.pred=b.1.event.bs, h0.truth=h0.1.truth.event)$estimate
  
  #1-transition-BP:
  b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
  b.1.bp=predict(b.1.event.bp,y1)
  # b.1.bp=BP.basis.finder(n,m1,y1)
  # n.event.1=length(y1[delta1==1])
  # b.1.bp.event=BP.basis.finder(n.event.1,m1,y1[delta1==1])
  phi.1.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m1+1)),spline.pred=b.1.event.bp,h0.truth=h0.1.truth.event)$estimate
  start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength, method="nls")
  
  # #2 transition:
  h0.2.truth.event <- alpha2.true*kappa2.true*y2[delta2==1]^(alpha2.true-1)
  ## 2-transition-BS:
  knot.loc.2 <- quantile(y2[delta2==1], ppoints(n.internalKnots.2))
  b.2.event.bs <- bSpline(y2[delta2==1], knots = knot.loc.2, degree = Bspline.degree.2, intercept=TRUE, Boundary.knots=bdy.knots.b.2)
  b.2.bs <- predict(b.2.event.bs, y2)
  phi.2.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.2), spline.pred=b.2.event.bs, h0.truth=h0.2.truth.event)$estimate
  #2-transition-BP:
  b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
  b.2.bp=predict(
    b.2.event.bp,y2)
  # b.2.bp=BP.basis.finder(n,m2,y2)
  # n.event.2=length(y2[delta2==1])
  # b.2.bp.event=BP.basis.finder(n.event.2,m1,y2[delta2==1])
  phi.2.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m2+1)),spline.pred=b.2.event.bp,h0.truth=h0.2.truth.event)$estimate
  start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength=500, method="nls")
  
  ## 3-transition:
  h0.3.truth.event <- alpha3.true*kappa3.true*(y2-y1)[delta1==1]^(alpha3.true-1)
  
  # #3 transition-BS:
  knot.loc.3 <- quantile((y2-y1)[delta1==1], ppoints(n.internalKnots.3))
  b.3.event.bs <- bSpline((y2-y1)[delta1==1], knots = knot.loc.3, degree = Bspline.degree.3, intercept=TRUE, Boundary.knots=bdy.knots.b.3.y2my1)
  b.3.y2my1.bs <- predict(b.3.event.bs, y2-y1)
  phi.3.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.3), spline.pred=b.3.event.bs, h0.truth=h0.3.truth.event)$estimate
  
  #3-transition-BP:
  # b.3.bp=BP.basis.finder(n,m3,(y2-y1))
  # n.event.3=length((y2-y1)[delta1==1])
  # b.3.bp.event=BP.basis.finder(n.event.3,m3,(y2-y1)[delta1==1])
  b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
  b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
  phi.3.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m3+1)),spline.pred=b.3.event.bp,h0.truth=h0.3.truth.event)$estimate
  start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength, method="nls")
  
  ## Start values for b-spline simulations    
  startVals.bs = c(phi.1.truth.bs,
                   phi.2.truth.bs,
                   phi.3.truth.bs)
  if (frailty == T) startVals.bs = c(startVals.bs, log(theta.true))
  startVals.bs = c(startVals.bs, beta1.true, beta2.true, beta3.true)
  
  
  beta1.start.bp=start1.bp[(length(phi.1.truth.bp)+1):length(start1.bp)]
  beta2.start.bp=start2.bp[(length(phi.2.truth.bp)+1):length(start2.bp)]
  beta3.start.bp=start3.bp[(length(phi.3.truth.bp)+1):length(start3.bp)]
  ## Start values for Bernstein simulations    
  startVals.bp = c(phi.1.truth.bp,
                   phi.2.truth.bp,
                   phi.3.truth.bp)
  if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true)+rnorm(1, mean=0, sd=0.1))
  # if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true))
  startVals.bp = c(startVals.bp, beta1.start.bp, beta2.start.bp, beta3.start.bp)
  
  
  phi.truth.BS=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  phi.truth.BP=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  
  #Weibull:
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  #BS:
  BSpara=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  if (frailty == T) BSpara=c(BSpara,log(theta.true))
  BSpara=c(BSpara,beta1.true,beta2.true,beta3.true)
  
  #BP:
  BPpara=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  if (frailty == T) BPpara=c(BPpara,log(theta.true))
  BPpara=c(BPpara,beta1.true,beta2.true,beta3.true)
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              startVals.bs = startVals.bs,
              startVals.bp=startVals.bp,
              frailty = frailty,
              b.1.bs=b.1.bs, b.2.bs=b.2.bs, b.3.y2my1.bs=b.3.y2my1.bs,
              bdy.knots.b.1=bdy.knots.b.1, 
              bdy.knots.b.2=bdy.knots.b.2, 
              bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1, 
              WBparamNames = WBparamNames,
              BSparamNames = BSparamNames,
              BPparamNames = BPparamNames,
              num.Bspline.params.1 = num.Bspline.params.1,
              num.Bspline.params.2 = num.Bspline.params.2,
              num.Bspline.params.3 = num.Bspline.params.3,
              m=m,
              WBpara = WBpara, 
              BSpara=BSpara,
              BPpara=BPpara,
              b.1.bp=b.1.bp,b.2.bp=b.2.bp,b.3.y2my1.bp=b.3.y2my1.bp,
              phi.truth.BS=phi.truth.BS,phi.truth.BP=phi.truth.BP))  
}  
genWB.simData.15cov.BP.BS<- function(weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution){
  #frailty parameter:
  # theta.true = ifelse(frailty == T, 0.25, 0)
  theta.true = ifelse(frailty == T, theta.given, 0)
  #Weibull parameters:
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
  
  
  
  m1=m[1]
  m2=m[2]
  m3=m[3]
  
  n.internalKnots.1 <- 1
  Bspline.degree.1 <- 1
  num.Bspline.params.1 <- (n.internalKnots.1) + (Bspline.degree.1 + 1)
  
  n.internalKnots.2 <- 1
  Bspline.degree.2 <- 1
  num.Bspline.params.2 <- (n.internalKnots.2) + (Bspline.degree.2 + 1)
  
  n.internalKnots.3 <- 1  ## Trying to add more flexibility - the previously estimated BH was off!
  Bspline.degree.3 <- 2
  num.Bspline.params.3 <- (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth))
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  
  
  num.Bspline.params.1 = (n.internalKnots.1) + (Bspline.degree.1 + 1)
  num.Bspline.params.2 = (n.internalKnots.2) + (Bspline.degree.2 + 1)
  num.Bspline.params.3 = (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  n.params.BS <- ifelse(frailty == T, num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true) + 1, 
                        num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true))
  
  BSparamNames <- c(paste0("eta1.", c(0:(m1))), paste0("eta2.", c(0:(m2))), paste0("eta3.", c(0:(m3))))
  BPparamNames=c(paste0("phi1.", c(0:(num.Bspline.params.1-1))), paste0("phi2.", c(0:(num.Bspline.params.2-1))), paste0("phi3.", c(0:(num.Bspline.params.3-1))))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    BSparamNames <- c(BSparamNames, "log(theta)")
    WBparamNames <- c(WBparamNames, "log(theta)") 
    BPparamNames=c(BPparamNames, "log(theta)")
  }
  BSparamNames <- c(BSparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  BPparamNames <- c(BPparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
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
  
  data0 = data.frame(first,second,third,fourth,fifth,sixth,seventh,eighth,ninth,tenth,eleventh,twelfth,thirteenth,fourteenth,fifteenth)
  
  
  
  # data0 = data.frame(first,second,third,fourth)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(lt1,lt2) # Origin starts at study entry
  
  Y.tmp = data.frame(simID.LT(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                              alpha1.true, alpha2.true, alpha3.true,
                              kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens,lt.type="unif", lt = lt,risks.status,frailty_distribution), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:20){
    Y[,i]=scale(Y[,i])
  }
  
  # Y.oracle=Y[,(1:10)]
  # data.oracle=data.frame(Y.oracle)
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  n=length(y1)
  bdy.knots.b.1 = c(0, max(y1)) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2-y1))
  
  # #1 tarnsition:
  h0.1.truth.event <- alpha1.true*kappa1.true*y1[delta1==1]^(alpha1.true-1)
  ## 1-transition-BS:
  knot.loc.1 <- quantile(y1[delta1==1], ppoints(n.internalKnots.1))
  b.1.event.bs <- bSpline(y1[delta1==1], knots = knot.loc.1, degree = Bspline.degree.1, intercept=TRUE, Boundary.knots=bdy.knots.b.1)
  b.1.bs <- predict(b.1.event.bs, y1)
  phi.1.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.1), spline.pred=b.1.event.bs, h0.truth=h0.1.truth.event)$estimate
  
  #1-transition-BP:
  b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
  b.1.bp=predict(b.1.event.bp,y1)
  # b.1.bp=BP.basis.finder(n,m1,y1)
  # n.event.1=length(y1[delta1==1])
  # b.1.bp.event=BP.basis.finder(n.event.1,m1,y1[delta1==1])
  phi.1.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m1+1)),spline.pred=b.1.event.bp,h0.truth=h0.1.truth.event)$estimate
  start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength, method="nls")
  
  # #2 transition:
  h0.2.truth.event <- alpha2.true*kappa2.true*y2[delta2==1]^(alpha2.true-1)
  ## 2-transition-BS:
  knot.loc.2 <- quantile(y2[delta2==1], ppoints(n.internalKnots.2))
  b.2.event.bs <- bSpline(y2[delta2==1], knots = knot.loc.2, degree = Bspline.degree.2, intercept=TRUE, Boundary.knots=bdy.knots.b.2)
  b.2.bs <- predict(b.2.event.bs, y2)
  phi.2.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.2), spline.pred=b.2.event.bs, h0.truth=h0.2.truth.event)$estimate
  #2-transition-BP:
  b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
  b.2.bp=predict(
    b.2.event.bp,y2)
  # b.2.bp=BP.basis.finder(n,m2,y2)
  # n.event.2=length(y2[delta2==1])
  # b.2.bp.event=BP.basis.finder(n.event.2,m1,y2[delta2==1])
  phi.2.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m2+1)),spline.pred=b.2.event.bp,h0.truth=h0.2.truth.event)$estimate
  start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength=500, method="nls")
  
  ## 3-transition:
  h0.3.truth.event <- alpha3.true*kappa3.true*(y2-y1)[delta1==1]^(alpha3.true-1)
  
  # #3 transition-BS:
  knot.loc.3 <- quantile((y2-y1)[delta1==1], ppoints(n.internalKnots.3))
  b.3.event.bs <- bSpline((y2-y1)[delta1==1], knots = knot.loc.3, degree = Bspline.degree.3, intercept=TRUE, Boundary.knots=bdy.knots.b.3.y2my1)
  b.3.y2my1.bs <- predict(b.3.event.bs, y2-y1)
  phi.3.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.3), spline.pred=b.3.event.bs, h0.truth=h0.3.truth.event)$estimate
  
  #3-transition-BP:
  # b.3.bp=BP.basis.finder(n,m3,(y2-y1))
  # n.event.3=length((y2-y1)[delta1==1])
  # b.3.bp.event=BP.basis.finder(n.event.3,m3,(y2-y1)[delta1==1])
  b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
  b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
  phi.3.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m3+1)),spline.pred=b.3.event.bp,h0.truth=h0.3.truth.event)$estimate
  start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength, method="nls")
  
  ## Start values for b-spline simulations    
  startVals.bs = c(phi.1.truth.bs,
                   phi.2.truth.bs,
                   phi.3.truth.bs)
  if (frailty == T) startVals.bs = c(startVals.bs, log(theta.true))
  startVals.bs = c(startVals.bs, beta1.true, beta2.true, beta3.true)
  
  
  beta1.start.bp=start1.bp[(length(phi.1.truth.bp)+1):length(start1.bp)]
  beta2.start.bp=start2.bp[(length(phi.2.truth.bp)+1):length(start2.bp)]
  beta3.start.bp=start3.bp[(length(phi.3.truth.bp)+1):length(start3.bp)]
  ## Start values for Bernstein simulations    
  startVals.bp = c(phi.1.truth.bp,
                   phi.2.truth.bp,
                   phi.3.truth.bp)
  if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true)+rnorm(1, mean=0, sd=0.1))
  # if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true))
  startVals.bp = c(startVals.bp, beta1.start.bp, beta2.start.bp, beta3.start.bp)
  
  
  phi.truth.BS=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  phi.truth.BP=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  
  #Weibull:
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  #BS:
  BSpara=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  if (frailty == T) BSpara=c(BSpara,log(theta.true))
  BSpara=c(BSpara,beta1.true,beta2.true,beta3.true)
  
  #BP:
  BPpara=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  if (frailty == T) BPpara=c(BPpara,log(theta.true))
  BPpara=c(BPpara,beta1.true,beta2.true,beta3.true)
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              startVals.bs = startVals.bs,
              startVals.bp=startVals.bp,
              frailty = frailty,
              b.1.bs=b.1.bs, b.2.bs=b.2.bs, b.3.y2my1.bs=b.3.y2my1.bs,
              bdy.knots.b.1=bdy.knots.b.1, 
              bdy.knots.b.2=bdy.knots.b.2, 
              bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1, 
              WBparamNames = WBparamNames,
              BSparamNames = BSparamNames,
              BPparamNames = BPparamNames,
              num.Bspline.params.1 = num.Bspline.params.1,
              num.Bspline.params.2 = num.Bspline.params.2,
              num.Bspline.params.3 = num.Bspline.params.3,
              m=m,
              WBpara = WBpara, 
              BSpara=BSpara,
              BPpara=BPpara,
              b.1.bp=b.1.bp,b.2.bp=b.2.bp,b.3.y2my1.bp=b.3.y2my1.bp,
              phi.truth.BS=phi.truth.BS,phi.truth.BP=phi.truth.BP))  
}  
genWB.simData.16cov.BP.BS<- function(weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution){
  #frailty parameter:
  # theta.true = ifelse(frailty == T, 0.25, 0)
  theta.true = ifelse(frailty == T, theta.given, 0)
  #Weibull parameters:
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
  
  
  
  m1=m[1]
  m2=m[2]
  m3=m[3]
  
  n.internalKnots.1 <- 1
  Bspline.degree.1 <- 1
  num.Bspline.params.1 <- (n.internalKnots.1) + (Bspline.degree.1 + 1)
  
  n.internalKnots.2 <- 1
  Bspline.degree.2 <- 1
  num.Bspline.params.2 <- (n.internalKnots.2) + (Bspline.degree.2 + 1)
  
  n.internalKnots.3 <- 1  ## Trying to add more flexibility - the previously estimated BH was off!
  Bspline.degree.3 <- 2
  num.Bspline.params.3 <- (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth))
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  
  
  num.Bspline.params.1 = (n.internalKnots.1) + (Bspline.degree.1 + 1)
  num.Bspline.params.2 = (n.internalKnots.2) + (Bspline.degree.2 + 1)
  num.Bspline.params.3 = (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  n.params.BS <- ifelse(frailty == T, num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true) + 1, 
                        num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true))
  
  BSparamNames <- c(paste0("eta1.", c(0:(m1))), paste0("eta2.", c(0:(m2))), paste0("eta3.", c(0:(m3))))
  BPparamNames=c(paste0("phi1.", c(0:(num.Bspline.params.1-1))), paste0("phi2.", c(0:(num.Bspline.params.2-1))), paste0("phi3.", c(0:(num.Bspline.params.3-1))))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    BSparamNames <- c(BSparamNames, "log(theta)")
    WBparamNames <- c(WBparamNames, "log(theta)") 
    BPparamNames=c(BPparamNames, "log(theta)")
  }
  BSparamNames <- c(BSparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  BPparamNames <- c(BPparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
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
  
  Y.tmp = data.frame(simID.LT(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                              alpha1.true, alpha2.true, alpha3.true,
                              kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens,lt.type="unif", lt = lt,risks.status,frailty_distribution), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:21){
    Y[,i]=scale(Y[,i])
  }
  
  # Y.oracle=Y[,(1:10)]
  # data.oracle=data.frame(Y.oracle)
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  n=length(y1)
  bdy.knots.b.1 = c(0, max(y1)) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2-y1))
  
  # #1 tarnsition:
  h0.1.truth.event <- alpha1.true*kappa1.true*y1[delta1==1]^(alpha1.true-1)
  ## 1-transition-BS:
  knot.loc.1 <- quantile(y1[delta1==1], ppoints(n.internalKnots.1))
  b.1.event.bs <- bSpline(y1[delta1==1], knots = knot.loc.1, degree = Bspline.degree.1, intercept=TRUE, Boundary.knots=bdy.knots.b.1)
  b.1.bs <- predict(b.1.event.bs, y1)
  phi.1.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.1), spline.pred=b.1.event.bs, h0.truth=h0.1.truth.event)$estimate
  
  #1-transition-BP:
  b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
  b.1.bp=predict(b.1.event.bp,y1)
  # b.1.bp=BP.basis.finder(n,m1,y1)
  # n.event.1=length(y1[delta1==1])
  # b.1.bp.event=BP.basis.finder(n.event.1,m1,y1[delta1==1])
  phi.1.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m1+1)),spline.pred=b.1.event.bp,h0.truth=h0.1.truth.event)$estimate
  start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength, method="nls")
  
  # #2 transition:
  h0.2.truth.event <- alpha2.true*kappa2.true*y2[delta2==1]^(alpha2.true-1)
  ## 2-transition-BS:
  knot.loc.2 <- quantile(y2[delta2==1], ppoints(n.internalKnots.2))
  b.2.event.bs <- bSpline(y2[delta2==1], knots = knot.loc.2, degree = Bspline.degree.2, intercept=TRUE, Boundary.knots=bdy.knots.b.2)
  b.2.bs <- predict(b.2.event.bs, y2)
  phi.2.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.2), spline.pred=b.2.event.bs, h0.truth=h0.2.truth.event)$estimate
  #2-transition-BP:
  b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
  b.2.bp=predict(
    b.2.event.bp,y2)
  # b.2.bp=BP.basis.finder(n,m2,y2)
  # n.event.2=length(y2[delta2==1])
  # b.2.bp.event=BP.basis.finder(n.event.2,m1,y2[delta2==1])
  phi.2.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m2+1)),spline.pred=b.2.event.bp,h0.truth=h0.2.truth.event)$estimate
  start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength=500, method="nls")
  
  ## 3-transition:
  h0.3.truth.event <- alpha3.true*kappa3.true*(y2-y1)[delta1==1]^(alpha3.true-1)
  
  # #3 transition-BS:
  knot.loc.3 <- quantile((y2-y1)[delta1==1], ppoints(n.internalKnots.3))
  b.3.event.bs <- bSpline((y2-y1)[delta1==1], knots = knot.loc.3, degree = Bspline.degree.3, intercept=TRUE, Boundary.knots=bdy.knots.b.3.y2my1)
  b.3.y2my1.bs <- predict(b.3.event.bs, y2-y1)
  phi.3.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.3), spline.pred=b.3.event.bs, h0.truth=h0.3.truth.event)$estimate
  
  #3-transition-BP:
  # b.3.bp=BP.basis.finder(n,m3,(y2-y1))
  # n.event.3=length((y2-y1)[delta1==1])
  # b.3.bp.event=BP.basis.finder(n.event.3,m3,(y2-y1)[delta1==1])
  b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
  b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
  phi.3.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m3+1)),spline.pred=b.3.event.bp,h0.truth=h0.3.truth.event)$estimate
  start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength, method="nls")
  
  ## Start values for b-spline simulations    
  startVals.bs = c(phi.1.truth.bs,
                   phi.2.truth.bs,
                   phi.3.truth.bs)
  if (frailty == T) startVals.bs = c(startVals.bs, log(theta.true))
  startVals.bs = c(startVals.bs, beta1.true, beta2.true, beta3.true)
  
  
  beta1.start.bp=start1.bp[(length(phi.1.truth.bp)+1):length(start1.bp)]
  beta2.start.bp=start2.bp[(length(phi.2.truth.bp)+1):length(start2.bp)]
  beta3.start.bp=start3.bp[(length(phi.3.truth.bp)+1):length(start3.bp)]
  ## Start values for Bernstein simulations    
  startVals.bp = c(phi.1.truth.bp,
                   phi.2.truth.bp,
                   phi.3.truth.bp)
  if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true)+rnorm(1, mean=0, sd=0.1))
  # if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true))
  startVals.bp = c(startVals.bp, beta1.start.bp, beta2.start.bp, beta3.start.bp)
  
  
  phi.truth.BS=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  phi.truth.BP=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  
  #Weibull:
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  #BS:
  BSpara=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  if (frailty == T) BSpara=c(BSpara,log(theta.true))
  BSpara=c(BSpara,beta1.true,beta2.true,beta3.true)
  
  #BP:
  BPpara=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  if (frailty == T) BPpara=c(BPpara,log(theta.true))
  BPpara=c(BPpara,beta1.true,beta2.true,beta3.true)
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              startVals.bs = startVals.bs,
              startVals.bp=startVals.bp,
              frailty = frailty,
              b.1.bs=b.1.bs, b.2.bs=b.2.bs, b.3.y2my1.bs=b.3.y2my1.bs,
              bdy.knots.b.1=bdy.knots.b.1, 
              bdy.knots.b.2=bdy.knots.b.2, 
              bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1, 
              WBparamNames = WBparamNames,
              BSparamNames = BSparamNames,
              BPparamNames = BPparamNames,
              num.Bspline.params.1 = num.Bspline.params.1,
              num.Bspline.params.2 = num.Bspline.params.2,
              num.Bspline.params.3 = num.Bspline.params.3,
              m=m,
              WBpara = WBpara, 
              BSpara=BSpara,
              BPpara=BPpara,
              b.1.bp=b.1.bp,b.2.bp=b.2.bp,b.3.y2my1.bp=b.3.y2my1.bp,
              phi.truth.BS=phi.truth.BS,phi.truth.BP=phi.truth.BP))  
}  

genWB.simData.20cov.BP.BS<- function(weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution){
  #frailty parameter:
  # theta.true = ifelse(frailty == T, 0.25, 0)
  theta.true = ifelse(frailty == T, theta.given, 0)
  #Weibull parameters:
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
  
  
  
  m1=m[1]
  m2=m[2]
  m3=m[3]
  
  n.internalKnots.1 <- 1
  Bspline.degree.1 <- 1
  num.Bspline.params.1 <- (n.internalKnots.1) + (Bspline.degree.1 + 1)
  
  n.internalKnots.2 <- 1
  Bspline.degree.2 <- 1
  num.Bspline.params.2 <- (n.internalKnots.2) + (Bspline.degree.2 + 1)
  
  n.internalKnots.3 <- 1  ## Trying to add more flexibility - the previously estimated BH was off!
  Bspline.degree.3 <- 2
  num.Bspline.params.3 <- (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth+twentyth), 
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth+twentyth),
                  as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth+thirteenth+fourteenth+fifteenth+sixteenth+seventeenth+eighteenth+ninteenth+twentyth))
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  
  
  num.Bspline.params.1 = (n.internalKnots.1) + (Bspline.degree.1 + 1)
  num.Bspline.params.2 = (n.internalKnots.2) + (Bspline.degree.2 + 1)
  num.Bspline.params.3 = (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  n.params.BS <- ifelse(frailty == T, num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true) + 1, 
                        num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true))
  
  BSparamNames <- c(paste0("eta1.", c(0:(m1))), paste0("eta2.", c(0:(m2))), paste0("eta3.", c(0:(m3))))
  BPparamNames=c(paste0("phi1.", c(0:(num.Bspline.params.1-1))), paste0("phi2.", c(0:(num.Bspline.params.2-1))), paste0("phi3.", c(0:(num.Bspline.params.3-1))))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    BSparamNames <- c(BSparamNames, "log(theta)")
    WBparamNames <- c(WBparamNames, "log(theta)") 
    BPparamNames=c(BPparamNames, "log(theta)")
  }
  BSparamNames <- c(BSparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  BPparamNames <- c(BPparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
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
  
  Y.tmp = data.frame(simID.LT(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                              alpha1.true, alpha2.true, alpha3.true,
                              kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens,lt.type="unif", lt = lt,risks.status,frailty_distribution), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:15){
    Y[,i]=scale(Y[,i])
  }
  
  # Y.oracle=Y[,(1:10)]
  # data.oracle=data.frame(Y.oracle)
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  n=length(y1)
  bdy.knots.b.1 = c(0, max(y1)) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2-y1))
  
  # #1 tarnsition:
  h0.1.truth.event <- alpha1.true*kappa1.true*y1[delta1==1]^(alpha1.true-1)
  ## 1-transition-BS:
  knot.loc.1 <- quantile(y1[delta1==1], ppoints(n.internalKnots.1))
  b.1.event.bs <- bSpline(y1[delta1==1], knots = knot.loc.1, degree = Bspline.degree.1, intercept=TRUE, Boundary.knots=bdy.knots.b.1)
  b.1.bs <- predict(b.1.event.bs, y1)
  phi.1.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.1), spline.pred=b.1.event.bs, h0.truth=h0.1.truth.event)$estimate
  
  #1-transition-BP:
  b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
  b.1.bp=predict(b.1.event.bp,y1)
  # b.1.bp=BP.basis.finder(n,m1,y1)
  # n.event.1=length(y1[delta1==1])
  # b.1.bp.event=BP.basis.finder(n.event.1,m1,y1[delta1==1])
  phi.1.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m1+1)),spline.pred=b.1.event.bp,h0.truth=h0.1.truth.event)$estimate
  start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength, method="nls")
  
  # #2 transition:
  h0.2.truth.event <- alpha2.true*kappa2.true*y2[delta2==1]^(alpha2.true-1)
  ## 2-transition-BS:
  knot.loc.2 <- quantile(y2[delta2==1], ppoints(n.internalKnots.2))
  b.2.event.bs <- bSpline(y2[delta2==1], knots = knot.loc.2, degree = Bspline.degree.2, intercept=TRUE, Boundary.knots=bdy.knots.b.2)
  b.2.bs <- predict(b.2.event.bs, y2)
  phi.2.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.2), spline.pred=b.2.event.bs, h0.truth=h0.2.truth.event)$estimate
  #2-transition-BP:
  b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
  b.2.bp=predict(
    b.2.event.bp,y2)
  # b.2.bp=BP.basis.finder(n,m2,y2)
  # n.event.2=length(y2[delta2==1])
  # b.2.bp.event=BP.basis.finder(n.event.2,m1,y2[delta2==1])
  phi.2.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m2+1)),spline.pred=b.2.event.bp,h0.truth=h0.2.truth.event)$estimate
  start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength=500, method="nls")
  
  ## 3-transition:
  h0.3.truth.event <- alpha3.true*kappa3.true*(y2-y1)[delta1==1]^(alpha3.true-1)
  
  # #3 transition-BS:
  knot.loc.3 <- quantile((y2-y1)[delta1==1], ppoints(n.internalKnots.3))
  b.3.event.bs <- bSpline((y2-y1)[delta1==1], knots = knot.loc.3, degree = Bspline.degree.3, intercept=TRUE, Boundary.knots=bdy.knots.b.3.y2my1)
  b.3.y2my1.bs <- predict(b.3.event.bs, y2-y1)
  phi.3.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.3), spline.pred=b.3.event.bs, h0.truth=h0.3.truth.event)$estimate
  
  #3-transition-BP:
  # b.3.bp=BP.basis.finder(n,m3,(y2-y1))
  # n.event.3=length((y2-y1)[delta1==1])
  # b.3.bp.event=BP.basis.finder(n.event.3,m3,(y2-y1)[delta1==1])
  b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
  b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
  phi.3.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m3+1)),spline.pred=b.3.event.bp,h0.truth=h0.3.truth.event)$estimate
  start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength, method="nls")
  
  ## Start values for b-spline simulations    
  startVals.bs = c(phi.1.truth.bs,
                   phi.2.truth.bs,
                   phi.3.truth.bs)
  if (frailty == T) startVals.bs = c(startVals.bs, log(theta.true))
  startVals.bs = c(startVals.bs, beta1.true, beta2.true, beta3.true)
  
  
  beta1.start.bp=start1.bp[(length(phi.1.truth.bp)+1):length(start1.bp)]
  beta2.start.bp=start2.bp[(length(phi.2.truth.bp)+1):length(start2.bp)]
  beta3.start.bp=start3.bp[(length(phi.3.truth.bp)+1):length(start3.bp)]
  ## Start values for Bernstein simulations    
  startVals.bp = c(phi.1.truth.bp,
                   phi.2.truth.bp,
                   phi.3.truth.bp)
  if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true)+rnorm(1, mean=0, sd=0.1))
  # if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true))
  startVals.bp = c(startVals.bp, beta1.start.bp, beta2.start.bp, beta3.start.bp)
  
  
  phi.truth.BS=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  phi.truth.BP=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  
  #Weibull:
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  #BS:
  BSpara=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  if (frailty == T) BSpara=c(BSpara,log(theta.true))
  BSpara=c(BSpara,beta1.true,beta2.true,beta3.true)
  
  #BP:
  BPpara=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  if (frailty == T) BPpara=c(BPpara,log(theta.true))
  BPpara=c(BPpara,beta1.true,beta2.true,beta3.true)
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              startVals.bs = startVals.bs,
              startVals.bp=startVals.bp,
              frailty = frailty,
              b.1.bs=b.1.bs, b.2.bs=b.2.bs, b.3.y2my1.bs=b.3.y2my1.bs,
              bdy.knots.b.1=bdy.knots.b.1, 
              bdy.knots.b.2=bdy.knots.b.2, 
              bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1, 
              WBparamNames = WBparamNames,
              BSparamNames = BSparamNames,
              BPparamNames = BPparamNames,
              num.Bspline.params.1 = num.Bspline.params.1,
              num.Bspline.params.2 = num.Bspline.params.2,
              num.Bspline.params.3 = num.Bspline.params.3,
              m=m,
              WBpara = WBpara, 
              BSpara=BSpara,
              BPpara=BPpara,
              b.1.bp=b.1.bp,b.2.bp=b.2.bp,b.3.y2my1.bp=b.3.y2my1.bp,
              phi.truth.BS=phi.truth.BS,phi.truth.BP=phi.truth.BP))  
}  

genWB.simData.10cov.BP.BS.grp<- function(G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n ,m,yvalslength, rho,frailty ,risks.status,theta.given,frailty_distribution){
  #frailty parameter:
  # theta.true = ifelse(frailty == T, 0.25, 0)
  theta.true = ifelse(frailty == T, theta.given, 0)
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
  
  
  cens = c(c1,c2)
  
  
  
  m1=m[1]
  m2=m[2]
  m3=m[3]
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  # #Trying:
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  ## B-spline baseline hazard specifications
  n.internalKnots.1 <- 1
  Bspline.degree.1 <- 1
  num.Bspline.params.1 <- (n.internalKnots.1) + (Bspline.degree.1 + 1)
  
  n.internalKnots.2 <- 1
  Bspline.degree.2 <- 1
  num.Bspline.params.2 <- (n.internalKnots.2) + (Bspline.degree.2 + 1)
  
  n.internalKnots.3 <- 1  ## Trying to add more flexibility - the previously estimated BH was off!
  Bspline.degree.3 <- 2
  num.Bspline.params.3 <- (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  
  lin.pred = list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth))
  
  if (risks.status=="fullyshared"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  }
  if (risks.status=="overlapping"){
    lin.pred.oracle = list(as.formula(~first+second+third+fourth), as.formula(~fourth+fifth+sixth+seventh), as.formula(~third+fourth+fifth+sixth))
  }
  # lin.pred = list(as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth), as.formula(~first+second+third+fourth))
  num.Bspline.params.1 = (n.internalKnots.1) + (Bspline.degree.1 + 1)
  num.Bspline.params.2 = (n.internalKnots.2) + (Bspline.degree.2 + 1)
  num.Bspline.params.3 = (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  n.params.BS <- ifelse(frailty == T, num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true) + 1, 
                        num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true))
  
  BSparamNames <- c(paste0("eta1.", c(0:(m1))), paste0("eta2.", c(0:(m2))), paste0("eta3.", c(0:(m3))))
  BPparamNames=c(paste0("phi1.", c(0:(num.Bspline.params.1-1))), paste0("phi2.", c(0:(num.Bspline.params.2-1))), paste0("phi3.", c(0:(num.Bspline.params.3-1))))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    BSparamNames <- c(BSparamNames, "log(theta)")
    WBparamNames <- c(WBparamNames, "log(theta)") 
    BPparamNames=c(BPparamNames, "log(theta)")
  }
  BSparamNames <- c(BSparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  BPparamNames <- c(BPparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
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
  
  lt = c(lt1,lt2) # Origin starts at study entry
  
  
  Y.tmp = data.frame(simID.LT(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                              alpha1.true, alpha2.true, alpha3.true,
                              kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens,lt.type="unif", lt = lt,risks.status,frailty_distribution), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  for (i in 6:15){
    Y[,i]=scale(Y[,i])
  }
  
  
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  n=length(y1)
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  
  
  
  
  n=length(y1)
  ## Estimate B-spline params
  ## log h0(t) = phi0 * B0(t) + ... + phik * Bk(t)
  bdy.knots.b.1 = c(0, max(y1)) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2-y1))
  
  # #1 tarnsition:
  h0.1.truth.event <- alpha1.true*kappa1.true*y1[delta1==1]^(alpha1.true-1)
  
  ## 1-transition-BS:
  knot.loc.1 <- quantile(y1[delta1==1], ppoints(n.internalKnots.1))
  b.1.event.bs <- bSpline(y1[delta1==1], knots = knot.loc.1, degree = Bspline.degree.1, intercept=TRUE, Boundary.knots=bdy.knots.b.1)
  b.1.bs <- predict(b.1.event.bs, y1)
  phi.1.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.1), spline.pred=b.1.event.bs, h0.truth=h0.1.truth.event)$estimate
  
  #1-transition-BP:
  b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
  b.1.bp=predict(b.1.event.bp,y1)
  # b.1.bp=BP.basis.finder(n,m1,y1)
  # n.event.1=length(y1[delta1==1])
  # b.1.bp.event=BP.basis.finder(n.event.1,m1,y1[delta1==1])
  phi.1.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m1+1)),spline.pred=b.1.event.bp,h0.truth=h0.1.truth.event)$estimate
  start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength, method="nls")
  
  # #2 transition:
  h0.2.truth.event <- alpha2.true*kappa2.true*y2[delta2==1]^(alpha2.true-1)
  
  ## 2-transition-BS:
  knot.loc.2 <- quantile(y2[delta2==1], ppoints(n.internalKnots.2))
  b.2.event.bs <- bSpline(y2[delta2==1], knots = knot.loc.2, degree = Bspline.degree.2, intercept=TRUE, Boundary.knots=bdy.knots.b.2)
  b.2.bs <- predict(b.2.event.bs, y2)
  phi.2.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.2), spline.pred=b.2.event.bs, h0.truth=h0.2.truth.event)$estimate
  
  #2-transition-BP:
  b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
  b.2.bp=predict(b.2.event.bp,y2)
  # b.2.bp=BP.basis.finder(n,m2,y2)
  # n.event.2=length(y2[delta2==1])
  # b.2.bp.event=BP.basis.finder(n.event.2,m1,y2[delta2==1])
  phi.2.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m2+1)),spline.pred=b.2.event.bp,h0.truth=h0.2.truth.event)$estimate
  start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength, method="nls")
  
  ## 3-transition:
  h0.3.truth.event <- alpha3.true*kappa3.true*(y2-y1)[delta1==1]^(alpha3.true-1)
  
  # #3 transition-BS:
  knot.loc.3 <- quantile((y2-y1)[delta1==1], ppoints(n.internalKnots.3))
  b.3.event.bs <- bSpline((y2-y1)[delta1==1], knots = knot.loc.3, degree = Bspline.degree.3, intercept=TRUE, Boundary.knots=bdy.knots.b.3.y2my1)
  b.3.y2my1.bs <- predict(b.3.event.bs, y2-y1)
  phi.3.truth.bs=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.3), spline.pred=b.3.event.bs, h0.truth=h0.3.truth.event)$estimate
  #3-transition-BP:
  # b.3.bp=BP.basis.finder(n,m3,(y2-y1))
  # n.event.3=length((y2-y1)[delta1==1])
  # b.3.bp.event=BP.basis.finder(n.event.3,m3,(y2-y1)[delta1==1])
  b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
  b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
  phi.3.truth.bp=nlm(f=obj.fn.1,p=rep(-1,(m3+1)),spline.pred=b.3.event.bp,h0.truth=h0.3.truth.event)$estimate
  start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength, method="nls")
  
  ## Start values for b-spline simulations    
  startVals.bs = c(phi.1.truth.bs,
                   phi.2.truth.bs,
                   phi.3.truth.bs)
  if (frailty == T) startVals.bs = c(startVals.bs, log(theta.true))
  startVals.bs = c(startVals.bs, beta1.true, beta2.true, beta3.true)
  
  beta1.start.bp=start1.bp[(length(phi.1.truth.bp)+1):length(start1.bp)]
  beta2.start.bp=start2.bp[(length(phi.2.truth.bp)+1):length(start2.bp)]
  beta3.start.bp=start3.bp[(length(phi.3.truth.bp)+1):length(start3.bp)]
  ## Start values for Bernstein simulations    
  startVals.bp = c(phi.1.truth.bp,
                   phi.2.truth.bp,
                   phi.3.truth.bp)
  if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true)+rnorm(1, mean=0, sd=0.1))
  # if (frailty == T) startVals.bp = c(startVals.bp, log(theta.true))
  startVals.bp = c(startVals.bp, beta1.start.bp, beta2.start.bp, beta3.start.bp)
  
  
  phi.truth.BS=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  phi.truth.BP=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  
  #Weibull:
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  #BS:
  BSpara=c(phi.1.truth.bs,phi.2.truth.bs,phi.3.truth.bs)
  if (frailty == T) BSpara=c(BSpara,log(theta.true))
  BSpara=c(BSpara,beta1.true,beta2.true,beta3.true)
  
  #BP:
  BPpara=c(phi.1.truth.bp,phi.2.truth.bp,phi.3.truth.bp)
  if (frailty == T) BPpara=c(BPpara,log(theta.true))
  BPpara=c(BPpara,beta1.true,beta2.true,beta3.true)
  
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              lin.pred.oracle=lin.pred.oracle,
              data = data,
              startVals.bs = startVals.bs,
              startVals.bp=startVals.bp,
              frailty = frailty,
              b.1.bs=b.1.bs, b.2.bs=b.2.bs, b.3.y2my1.bs=b.3.y2my1.bs,
              bdy.knots.b.1=bdy.knots.b.1, 
              bdy.knots.b.2=bdy.knots.b.2, 
              bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1, 
              WBparamNames = WBparamNames,
              BSparamNames = BSparamNames,
              BPparamNames = BPparamNames,
              num.Bspline.params.1 = num.Bspline.params.1,
              num.Bspline.params.2 = num.Bspline.params.2,
              num.Bspline.params.3 = num.Bspline.params.3,
              m=m,
              WBpara = WBpara, 
              BSpara=BSpara,
              BPpara=BPpara,
              b.1.bp=b.1.bp,b.2.bp=b.2.bp,b.3.y2my1.bp=b.3.y2my1.bp,
              phi.truth.BS=phi.truth.BS,phi.truth.BP=phi.truth.BP)) 
}  


# Bernstein poly, and B-splines ----------------------------------------------------------------


logLike.SCR.SM.LT.bSpline.bp.dropPrevCases <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE,
                                                       b.1,  
                                                       b.2,  
                                                       b.3.y2my1)
{
  ##
  if (is.vector(Xmat1)==T) Xmat1 = matrix(Xmat1, nrow=1)
  if (is.vector(Xmat2)==T) Xmat2 = matrix(Xmat2, nrow=1)
  if (is.vector(Xmat3)==T) Xmat3 = matrix(Xmat3, nrow=1)
  num.Bspline.params.1 <- ncol(b.1)
  num.Bspline.params.2 <- ncol(b.2)
  num.Bspline.params.3 <- ncol(b.3.y2my1)
  phi1 <- para[1:(1+num.Bspline.params.1-1)]
  phi2 <- para[(1+num.Bspline.params.1):(1+num.Bspline.params.1 + num.Bspline.params.2 - 1)]
  phi3 <- para[(1+num.Bspline.params.1 + num.Bspline.params.2):(1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)]
  
  if(frailty == TRUE){
    theta    <- exp(para[(2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, (2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1), (1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1))
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
  ##
  ## log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  ## Uses h1 - need to bSpline
  ## log.h1,0(t) = phi0 * B0(t) + ... + phik * Bk(t)
  B0.1.y1 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(b.1)))
  log.h1star.y1 <- B0.1.y1 + eta.1
  
  ## log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  B0.2.y1 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, y1))))
  log.h2star.y1 <- B0.2.y1 + eta.2
  
  ## log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  B0.2.y2 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(b.2)))
  log.h2star.y2 <- B0.2.y2 + eta.2
  
  ## log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  B0.3.y2my1 <- as.vector(matrix(phi3, nrow=1) %*% t(as.matrix(b.3.y2my1)))
  log.h3star.y2 <- B0.3.y2my1 + eta.3
  
  ## q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  ## Uses h1 - needs bSpline
  ## Essentially uses H0.1.y1 and H0.2.y1
  ## exp(phi.3.truth %*% t(b.3.event))
  
  h0.1.y1 <- exp(B0.1.y1)
  h0.1.y1.interpolate.func <- approxfun(y1, h0.1.y1, rule=2)
  H0.1.y1 <- sapply(y1, function(x) integrate(h0.1.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y1 <- exp(B0.2.y1)
  h0.2.y1.interpolate.func <- approxfun(y1, h0.2.y1, rule=2)
  H0.2.y1 <- sapply(y1, function(x) integrate(h0.2.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y1 <- H0.1.y1 * exp(eta.1) + H0.2.y1 * exp(eta.2)
  
  ## q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.y2 and H0.2.y2
  
  B0.1.y2 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, y2))))
  h0.1.y2 <- exp(B0.1.y2)
  h0.1.y2.interpolate.func <- approxfun(y2, h0.1.y2, rule=2)
  H0.1.y2 <- sapply(y2, function(x) integrate(h0.1.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y2 <- exp(B0.2.y2)
  h0.2.y2.interpolate.func <- approxfun(y2, h0.2.y2, rule=2)
  H0.2.y2 <- sapply(y2, function(x) integrate(h0.2.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y2 <- H0.1.y2 * exp(eta.1) + H0.2.y2 * exp(eta.2)
  
  ## q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.l and H0.2.l
  
  B0.1.l <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, l))))
  h0.1.l <- exp(B0.1.l)
  h0.1.l.interpolate.func <- approxfun(l, h0.1.l, rule=2)
  H0.1.l <- sapply(l, function(x) integrate(h0.1.l.interpolate.func, 0, x, stop.on.error = F)$value)
  
  B0.2.l <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, l))))
  h0.2.l <- exp(B0.2.l)
  h0.2.l.interpolate.func <- approxfun(l, h0.2.l, rule=2)
  H0.2.l <- sapply(l, function(x) integrate(h0.2.l.interpolate.func, 0, x, stop.on.error = F)$value)
  q.l <- H0.1.l * exp(eta.1) + H0.2.l * exp(eta.2)
  
  ## w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  ## This is H0.3(y2-y1)
  
  h0.3.y2my1 <- exp(B0.3.y2my1)
  h0.3.y2my1.interpolate.func <- approxfun(y2-y1, h0.3.y2my1, rule=2)
  H0.3.y2my1 <- sapply(y2-y1, function(x) integrate(h0.3.y2my1.interpolate.func, 0, x, stop.on.error = F)$value)
  w.y1.y2 <- H0.3.y2my1 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (k1 + k2.y1))))
    logLike2 <- log.h2star.y1 - ((thetaInv + 1) * log(1 + (theta * k2.y1)))  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (k1 + k2.y1))))
    logLike4 <- - thetaInv * log(1 + (theta * k2.y1))  ## Making in terms of y1
  }
  if(frailty == FALSE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
    logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - (k1 + k2.y1)
    logLike4 <- - k2.y1  ## Making in terms of y1
  }
  ##
  loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) 
  return(-loglh)
}


logLike.SCR.SM.LT.bSpline.bp.dropPrevCases.beta.fixed <- function(para,beta1.from.varsel,beta2.from.varsel,beta3.from.varsel, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE,
                                                       b.1,  
                                                       b.2,  
                                                       b.3.y2my1)
{
  ##
  if (is.vector(Xmat1)==T) Xmat1 = matrix(Xmat1, nrow=1)
  if (is.vector(Xmat2)==T) Xmat2 = matrix(Xmat2, nrow=1)
  if (is.vector(Xmat3)==T) Xmat3 = matrix(Xmat3, nrow=1)
  num.Bspline.params.1 <- ncol(b.1)
  num.Bspline.params.2 <- ncol(b.2)
  num.Bspline.params.3 <- ncol(b.3.y2my1)
  phi1 <- para[1:(1+num.Bspline.params.1-1)]
  phi2 <- para[(1+num.Bspline.params.1):(1+num.Bspline.params.1 + num.Bspline.params.2 - 1)]
  phi3 <- para[(1+num.Bspline.params.1 + num.Bspline.params.2):(1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)]
  
  if(frailty == TRUE){
    theta    <- exp(para[(2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, (2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1), (1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1))
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% beta1.from.varsel)
  eta.2 <- as.vector(Xmat2 %*% beta2.from.varsel)
  eta.3 <- as.vector(Xmat3 %*% beta3.from.varsel)
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1 & l < y1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1 & l < y1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0 & l < y1)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0 & l < y1)
  ##
  ## log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  ## Uses h1 - need to bSpline
  ## log.h1,0(t) = phi0 * B0(t) + ... + phik * Bk(t)
  B0.1.y1 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(b.1)))
  log.h1star.y1 <- B0.1.y1 + eta.1
  
  ## log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  B0.2.y1 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, y1))))
  log.h2star.y1 <- B0.2.y1 + eta.2
  
  ## log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  B0.2.y2 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(b.2)))
  log.h2star.y2 <- B0.2.y2 + eta.2
  
  ## log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  B0.3.y2my1 <- as.vector(matrix(phi3, nrow=1) %*% t(as.matrix(b.3.y2my1)))
  log.h3star.y2 <- B0.3.y2my1 + eta.3
  
  ## q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  ## Uses h1 - needs bSpline
  ## Essentially uses H0.1.y1 and H0.2.y1
  ## exp(phi.3.truth %*% t(b.3.event))
  
  h0.1.y1 <- exp(B0.1.y1)
  h0.1.y1.interpolate.func <- approxfun(y1, h0.1.y1, rule=2)
  H0.1.y1 <- sapply(y1, function(x) integrate(h0.1.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y1 <- exp(B0.2.y1)
  h0.2.y1.interpolate.func <- approxfun(y1, h0.2.y1, rule=2)
  H0.2.y1 <- sapply(y1, function(x) integrate(h0.2.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y1 <- H0.1.y1 * exp(eta.1) + H0.2.y1 * exp(eta.2)
  
  ## q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.y2 and H0.2.y2
  
  B0.1.y2 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, y2))))
  h0.1.y2 <- exp(B0.1.y2)
  h0.1.y2.interpolate.func <- approxfun(y2, h0.1.y2, rule=2)
  H0.1.y2 <- sapply(y2, function(x) integrate(h0.1.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y2 <- exp(B0.2.y2)
  h0.2.y2.interpolate.func <- approxfun(y2, h0.2.y2, rule=2)
  H0.2.y2 <- sapply(y2, function(x) integrate(h0.2.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y2 <- H0.1.y2 * exp(eta.1) + H0.2.y2 * exp(eta.2)
  
  ## q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.l and H0.2.l
  
  B0.1.l <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, l))))
  h0.1.l <- exp(B0.1.l)
  h0.1.l.interpolate.func <- approxfun(l, h0.1.l, rule=2)
  H0.1.l <- sapply(l, function(x) integrate(h0.1.l.interpolate.func, 0, x, stop.on.error = F)$value)
  
  B0.2.l <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, l))))
  h0.2.l <- exp(B0.2.l)
  h0.2.l.interpolate.func <- approxfun(l, h0.2.l, rule=2)
  H0.2.l <- sapply(l, function(x) integrate(h0.2.l.interpolate.func, 0, x, stop.on.error = F)$value)
  q.l <- H0.1.l * exp(eta.1) + H0.2.l * exp(eta.2)
  
  ## w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  ## This is H0.3(y2-y1)
  
  h0.3.y2my1 <- exp(B0.3.y2my1)
  h0.3.y2my1.interpolate.func <- approxfun(y2-y1, h0.3.y2my1, rule=2)
  H0.3.y2my1 <- sapply(y2-y1, function(x) integrate(h0.3.y2my1.interpolate.func, 0, x, stop.on.error = F)$value)
  w.y1.y2 <- H0.3.y2my1 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (k1 + k2.y1))))
    logLike2 <- log.h2star.y1 - ((thetaInv + 1) * log(1 + (theta * k2.y1)))  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (k1 + k2.y1))))
    logLike4 <- - thetaInv * log(1 + (theta * k2.y1))  ## Making in terms of y1
  }
  if(frailty == FALSE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
    logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - (k1 + k2.y1)
    logLike4 <- - k2.y1  ## Making in terms of y1
  }
  ##
  loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) 
  return(-loglh)
}
dlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE,
                                                        b.1,  
                                                        b.2,  
                                                        b.3.y2my1)
{
  ##
  if (is.vector(Xmat1)==T) Xmat1 = matrix(Xmat1, nrow=1)
  if (is.vector(Xmat2)==T) Xmat2 = matrix(Xmat2, nrow=1)
  if (is.vector(Xmat3)==T) Xmat3 = matrix(Xmat3, nrow=1)
  num.Bspline.params.1 <- ncol(b.1)
  num.Bspline.params.2 <- ncol(b.2)
  num.Bspline.params.3 <- ncol(b.3.y2my1)
  phi1 <- para[1:(1+num.Bspline.params.1-1)]
  phi2 <- para[(1+num.Bspline.params.1):(1+num.Bspline.params.1 + num.Bspline.params.2 - 1)]
  phi3 <- para[(1+num.Bspline.params.1 + num.Bspline.params.2):(1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)]
  
  if(frailty == TRUE){
    theta    <- exp(para[(2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, (2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1), (1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1))
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
  ##
  ## log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  ## Uses h1 - need to bSpline
  ## log.h1,0(t) = phi0 * B0(t) + ... + phik * Bk(t)
  B0.1.y1 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(b.1)))
  log.h1star.y1 <- B0.1.y1 + eta.1
  
  ## log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  B0.2.y1 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, y1))))
  log.h2star.y1 <- B0.2.y1 + eta.2
  
  ## log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  B0.2.y2 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(b.2)))
  log.h2star.y2 <- B0.2.y2 + eta.2
  
  ## log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  B0.3.y2my1 <- as.vector(matrix(phi3, nrow=1) %*% t(as.matrix(b.3.y2my1)))
  log.h3star.y2 <- B0.3.y2my1 + eta.3
  
  ## q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  ## Uses h1 - needs bSpline
  ## Essentially uses H0.1.y1 and H0.2.y1
  ## exp(phi.3.truth %*% t(b.3.event))
  
  h0.1.y1 <- exp(B0.1.y1)
  h0.1.y1.interpolate.func <- approxfun(y1, h0.1.y1, rule=2)
  H0.1.y1 <- sapply(y1, function(x) integrate(h0.1.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y1 <- exp(B0.2.y1)
  h0.2.y1.interpolate.func <- approxfun(y1, h0.2.y1, rule=2)
  H0.2.y1 <- sapply(y1, function(x) integrate(h0.2.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y1 <- H0.1.y1 * exp(eta.1) + H0.2.y1 * exp(eta.2)
  
  ## q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.y2 and H0.2.y2
  
  B0.1.y2 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, y2))))
  h0.1.y2 <- exp(B0.1.y2)
  h0.1.y2.interpolate.func <- approxfun(y2, h0.1.y2, rule=2)
  H0.1.y2 <- sapply(y2, function(x) integrate(h0.1.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y2 <- exp(B0.2.y2)
  h0.2.y2.interpolate.func <- approxfun(y2, h0.2.y2, rule=2)
  H0.2.y2 <- sapply(y2, function(x) integrate(h0.2.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y2 <- H0.1.y2 * exp(eta.1) + H0.2.y2 * exp(eta.2)
  
  ## q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.l and H0.2.l
  
  B0.1.l <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, l))))
  h0.1.l <- exp(B0.1.l)
  h0.1.l.interpolate.func <- approxfun(l, h0.1.l, rule=2)
  H0.1.l <- sapply(l, function(x) integrate(h0.1.l.interpolate.func, 0, x, stop.on.error = F)$value)
  
  B0.2.l <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, l))))
  h0.2.l <- exp(B0.2.l)
  h0.2.l.interpolate.func <- approxfun(l, h0.2.l, rule=2)
  H0.2.l <- sapply(l, function(x) integrate(h0.2.l.interpolate.func, 0, x, stop.on.error = F)$value)
  q.l <- H0.1.l * exp(eta.1) + H0.2.l * exp(eta.2)
  
  ## w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  ## This is H0.3(y2-y1)
  
  h0.3.y2my1 <- exp(B0.3.y2my1)
  h0.3.y2my1.interpolate.func <- approxfun(y2-y1, h0.3.y2my1, rule=2)
  H0.3.y2my1 <- sapply(y2-y1, function(x) integrate(h0.3.y2my1.interpolate.func, 0, x, stop.on.error = F)$value)
  w.y1.y2 <- H0.3.y2my1 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (k1 + k2.y1))))
    logLike2 <- log.h2star.y1 - ((thetaInv + 1) * log(1 + (theta * k2.y1)))  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (k1 + k2.y1))))
    logLike4 <- - thetaInv * log(1 + (theta * k2.y1))  ## Making in terms of y1
    n=NROW(eta.1)
    p=ncol(Xmat1)
    #score_ij-->i(l_i),(j:beta_j)--->score23:dl_2(\phi)/d\beta_3.
    #we have l_i; i=1,2,3,4, and beta_j; j=1,2,3,4.
    A1=(H0.1.y1*exp(eta.1))-(H0.1.l*exp(eta.1))
    denom1=1+(theta*(k1+k2.y1))
    denom2=1+(theta*k2.y1)
    A2=(H0.2.y1 *exp(eta.2))-(H0.2.l *exp(eta.2))
    A3=(H0.3.y2my1 *exp(eta.3))
    
    
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
  # if(frailty == FALSE)
  # {
  #   logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
  #   logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
  #   logLike3 <- log.h1star.y1 - (k1 + k2.y1)
  #   logLike4 <- - k2.y1  ## Making in terms of y1
  # }
  ##
  return(-score)
}
ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE,
                                                         b.1,  
                                                         b.2,  
                                                         b.3.y2my1)
{
  ##
  if (is.vector(Xmat1)==T) Xmat1 = matrix(Xmat1, nrow=1)
  if (is.vector(Xmat2)==T) Xmat2 = matrix(Xmat2, nrow=1)
  if (is.vector(Xmat3)==T) Xmat3 = matrix(Xmat3, nrow=1)
  num.Bspline.params.1 <- ncol(b.1)
  num.Bspline.params.2 <- ncol(b.2)
  num.Bspline.params.3 <- ncol(b.3.y2my1)
  phi1 <- para[1:(1+num.Bspline.params.1-1)]
  phi2 <- para[(1+num.Bspline.params.1):(1+num.Bspline.params.1 + num.Bspline.params.2 - 1)]
  phi3 <- para[(1+num.Bspline.params.1 + num.Bspline.params.2):(1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)]
  
  if(frailty == TRUE){
    theta    <- exp(para[(2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, (2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1), (1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1))
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
  ##
  ## log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  ## Uses h1 - need to bSpline
  ## log.h1,0(t) = phi0 * B0(t) + ... + phik * Bk(t)
  B0.1.y1 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(b.1)))
  log.h1star.y1 <- B0.1.y1 + eta.1
  
  ## log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  B0.2.y1 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, y1))))
  log.h2star.y1 <- B0.2.y1 + eta.2
  
  ## log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  B0.2.y2 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(b.2)))
  log.h2star.y2 <- B0.2.y2 + eta.2
  
  ## log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  B0.3.y2my1 <- as.vector(matrix(phi3, nrow=1) %*% t(as.matrix(b.3.y2my1)))
  log.h3star.y2 <- B0.3.y2my1 + eta.3
  
  ## q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  ## Uses h1 - needs bSpline
  ## Essentially uses H0.1.y1 and H0.2.y1
  ## exp(phi.3.truth %*% t(b.3.event))
  
  h0.1.y1 <- exp(B0.1.y1)
  h0.1.y1.interpolate.func <- approxfun(y1, h0.1.y1, rule=2)
  H0.1.y1 <- sapply(y1, function(x) integrate(h0.1.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y1 <- exp(B0.2.y1)
  h0.2.y1.interpolate.func <- approxfun(y1, h0.2.y1, rule=2)
  H0.2.y1 <- sapply(y1, function(x) integrate(h0.2.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y1 <- H0.1.y1 * exp(eta.1) + H0.2.y1 * exp(eta.2)
  
  ## q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.y2 and H0.2.y2
  
  B0.1.y2 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, y2))))
  h0.1.y2 <- exp(B0.1.y2)
  h0.1.y2.interpolate.func <- approxfun(y2, h0.1.y2, rule=2)
  H0.1.y2 <- sapply(y2, function(x) integrate(h0.1.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y2 <- exp(B0.2.y2)
  h0.2.y2.interpolate.func <- approxfun(y2, h0.2.y2, rule=2)
  H0.2.y2 <- sapply(y2, function(x) integrate(h0.2.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y2 <- H0.1.y2 * exp(eta.1) + H0.2.y2 * exp(eta.2)
  
  ## q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.l and H0.2.l
  
  B0.1.l <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, l))))
  h0.1.l <- exp(B0.1.l)
  h0.1.l.interpolate.func <- approxfun(l, h0.1.l, rule=2)
  H0.1.l <- sapply(l, function(x) integrate(h0.1.l.interpolate.func, 0, x, stop.on.error = F)$value)
  
  B0.2.l <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, l))))
  h0.2.l <- exp(B0.2.l)
  h0.2.l.interpolate.func <- approxfun(l, h0.2.l, rule=2)
  H0.2.l <- sapply(l, function(x) integrate(h0.2.l.interpolate.func, 0, x, stop.on.error = F)$value)
  q.l <- H0.1.l * exp(eta.1) + H0.2.l * exp(eta.2)
  
  ## w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  ## This is H0.3(y2-y1)
  
  h0.3.y2my1 <- exp(B0.3.y2my1)
  h0.3.y2my1.interpolate.func <- approxfun(y2-y1, h0.3.y2my1, rule=2)
  H0.3.y2my1 <- sapply(y2-y1, function(x) integrate(h0.3.y2my1.interpolate.func, 0, x, stop.on.error = F)$value)
  w.y1.y2 <- H0.3.y2my1 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  ##
  if(frailty == TRUE)
  {
    n=NROW(eta.1)
    p=ncol(Xmat1)
    #score_ij-->i(l_i),(j:beta_j)--->score23:dl_2(\phi)/d\beta_3.
    #we have l_i; i=1,2,3,4, and beta_j; j=1,2,3,4.
    A1=(theta*H0.1.y1*exp(eta.1))-(theta*H0.1.l *exp(eta.1))
    A2=(theta*H0.2.y1 *exp(eta.2))-(theta*H0.2.l *exp(eta.2))
    A3=theta*(H0.3.y2my1 *exp(eta.3))
    A4=(H0.1.y1-H0.1.l )*exp(eta.1)
    A5=(H0.1.y1*exp(eta.1))-(H0.1.l*exp(eta.1))
    A6=(H0.2.y1 *exp(eta.2))-(H0.2.l *exp(eta.2))
    A7=(H0.3.y2my1 *exp(eta.3))
    
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
  # if(frailty == FALSE)
  # {
  #   logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
  #   logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
  #   logLike3 <- log.h1star.y1 - (k1 + k2.y1)
  #   logLike4 <- - k2.y1  ## Making in terms of y1
  # }
  ##
  dscore=rbind(cbind(dscore11,dscore12,dscore13),cbind(dscore21,dscore22,dscore23),cbind(dscore31,dscore32,dscore33))
  ##
  return(-dscore)
}

FreqID.LT.bSpline.bp <- function(Y,  lin.pred, data, startVals, frailty, 
                                 b.1, b.2, b.3.y2my1, 
                                 bdy.knots.b.1, 
                                 bdy.knots.b.2, 
                                 bdy.knots.b.3.y2my1,
                                 method)
{	
  
  ##
  y1     <- as.vector(Y[,1])
  delta1 <- as.vector(Y[,2])
  y2     <- as.vector(Y[,3])
  delta2 <- as.vector(Y[,4])
  l      <- as.vector(Y[,5])
  Xmat1  <- as.matrix(model.frame(lin.pred[[1]], data=data))  
  Xmat2  <- as.matrix(model.frame(lin.pred[[2]], data=data))  
  Xmat3  <- as.matrix(model.frame(lin.pred[[3]], data=data))
  ##
  cat("This is going to take quite awhile...~30 minutes hour for ~300 observations")
  if (method == "optim"){
    logLike <- function(p) logLike.SCR.SM.LT.bSpline.bp.dropPrevCases(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
                                                                      b.1 = b.1,  
                                                                      b.2 = b.2, 
                                                                      b.3.y2my1 = b.3.y2my1)
    
    optim.control = list(REPORT = 50, maxit = 2000)  
    
    fit1 <- optim(startVals,   
                  logLike, hessian = TRUE, method="Nelder-Mead", control = optim.control)
    H = fit1$hessian
    est = fit1$par
    logLike=-fit1$value
    code = fit1$convergence}
  if (method == "nlm"){  
    
    fit2 <- nlm(logLike.SCR.SM.LT.bSpline.bp.dropPrevCases, p=startVals, y1=y1, y2=y2,
                delta1=delta1, delta2=delta2, l=l,
                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
                b.1 = b.1,
                b.2 = b.2,
                b.3.y2my1 = b.3.y2my1, hessian=T)
    H = fit2$hessian
    logLike = -fit2$minimum
    est = fit2$estimate
    code = fit2$code}
  ##
  myLabels <- c("phi1", "phi2", "phi3")
  if(frailty == TRUE) myLabels <- c(myLabels, "log(theta)")
  myLabels <- c(myLabels, colnames(Xmat1), colnames(Xmat2), colnames(Xmat3))
  nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
  
  value <- list(startingValues=startVals,estimate=est, H=H, logLike=logLike, myLabels=myLabels, frailty=frailty, nP=nP, 
                code = code)
  
  class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
  
  return(value)
  ##
  invisible()
}

FreqID.LT.bSpline.bp.oracle.overlapping <- function(Y, lin.pred.oracle, data.oracle.1,data.oracle.2,data.oracle.3,
                                                    startVals, frailty=TRUE, 
                                                    b.1, b.2, b.3.y2my1, 
                                                    bdy.knots.b.1, 
                                                    bdy.knots.b.2, 
                                                    bdy.knots.b.3.y2my1,
                                                    method )
{	
  
  ##
  y1     <- as.vector(Y[,1])
  delta1 <- as.vector(Y[,2])
  y2     <- as.vector(Y[,3])
  delta2 <- as.vector(Y[,4])
  l      <- as.vector(Y[,5])
  
  
  Xmat1  <- as.matrix(model.frame(lin.pred.oracle[[1]], data=data.oracle.1))
  Xmat2  <- as.matrix(model.frame(lin.pred.oracle[[2]], data=data.oracle.2))
  Xmat3  <- as.matrix(model.frame(lin.pred.oracle[[3]], data=data.oracle.3))
  
  ##
  cat("This is going to take quite awhile...~30 minutes hour for ~300 observations")
  if (method == "optim"){
    logLike <- function(p) logLike.SCR.SM.LT.bSpline.bp.dropPrevCases(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
                                                                      b.1 = b.1,  
                                                                      b.2 = b.2, 
                                                                      b.3.y2my1 = b.3.y2my1)
    
    optim.control = list(REPORT = 50, maxit = 2000)  
    
    fit1 <- optim(startVals,   
                  logLike, hessian = TRUE, method="Nelder-Mead", control = optim.control)
    H = fit1$hessian
    est = fit1$par
    logLike=-fit1$value
    code = fit1$convergence
  }
  if (method == "nlm"){ 
    
    fit2 <- nlm(logLike.SCR.SM.LT.bSpline.bp.dropPrevCases, p=startVals, y1=y1, y2=y2,
                delta1=delta1, delta2=delta2, l=l,
                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
                b.1 = b.1,
                b.2 = b.2,
                b.3.y2my1 = b.3.y2my1, hessian=T)
    H = fit2$hessian
    logLike = -fit2$minimum
    est = fit2$estimate
    code = fit2$code}
  ##
  myLabels <- c("phi1", "phi2", "phi3")
  if(frailty == TRUE) myLabels <- c(myLabels, "log(theta)")
  myLabels <- c(myLabels, colnames(Xmat1), colnames(Xmat2), colnames(Xmat3))
  nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
  
  value <- list(estimate=est, H=H, logLike=logLike, myLabels=myLabels, frailty=frailty, nP=nP, 
                code = code)
  
  class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
  
  return(value)
  ##
  invisible()
}


get.Bspline.startVals <- function(y, delta, lin.pred, data, b, knot.loc, Bspline.degree, method="nls"){
  ## 
  num.Bspline.params <- length(knot.loc) + Bspline.degree + 1
  
  fitCox <- coxph(as.formula(paste("Surv(y, delta)", paste(lin.pred, collapse = ""))), data=data)
  H0.info <- basehaz(fitCox, centered = F) ## Get cumulative hazard from Cox
  H0.est.linInterpol <- approxfun(H0.info$time, H0.info$hazard, rule = 2) ## Get lin interpolated function
  yvals <- seq(min(y), max(y), length=500) ## get plenty of y-values
  h0.est.func2 <- numdiff(H0.est.linInterpol, yvals) ## estimate h0 using numerical differentiation
  h0.est.smooth.xy <- loess.smooth(yvals, h0.est.func2, evaluation=500) ## Smooth out using loess
  
  ## Now using least squares to get estimated phi coefficients based on loess-smoothed h0
  h0.est.smooth <- h0.est.smooth.xy$y
  yvals.bspline <- data.frame(predict(b, yvals))
  
  obj.fn.1 <- function(eta.vec, spline.pred, h0.truth){ 
    ##
    spline.vals <- sapply(1:nrow(spline.pred), function(i) sum(spline.pred[i,]*(eta.vec)))
    sq.diffs <- (log(h0.truth)-spline.vals)^2
    return(sum(sq.diffs))
  }
  if (method == "nlm"){
    phi.h0.est.smooth=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params), 
                          spline.pred=yvals.bspline, 
                          h0.truth=h0.est.smooth)$estimate
    return(c(phi.h0.est.smooth, fitCox$coef))
  }
  if (method == "nls"){
    ## Alternatively -- this seems to give better results? 
    g <- function(bSpline.mat, phi.vec){
      phi.list <- as.list(phi.vec)
      mat.list <- as.list(data.frame(bSpline.mat))
      return(apply(do.call(rbind, Map('*', phi.list, mat.list)), 2, sum))
    }
    phi.h0.est.smooth.ls <- nls(h0.est.smooth ~ exp(g(yvals.bspline, phi.vec)), 
                                start = list(phi.vec = rep(-1, ncol(yvals.bspline))))
    phi.h0.est.smooth <- summary(phi.h0.est.smooth.ls)$parameters[,1]  ## These will be good for starting values!
    return(c(phi.h0.est.smooth, fitCox$coef))
  }
}

get.BP.startVals <- function(BPdegree,y, delta, lin.pred, data, b,yvalslength, method="nls"){
  ## 
  num.BP.params <- BPdegree+1
  
  fitCox <- coxph(as.formula(paste("Surv(y, delta)", paste(lin.pred, collapse = ""))), data=data)
  H0.info <- basehaz(fitCox, centered = F) ## Get cumulative hazard from Cox
  H0.est.linInterpol <- approxfun(H0.info$time, H0.info$hazard, rule = 2) ## Get lin interpolated function
  yvals <- seq(min(y), max(y), length=yvalslength) ## get plenty of y-values
  h0.est.func2 <- numdiff(H0.est.linInterpol, yvals) ## estimate h0 using numerical differentiation
  h0.est.smooth.xy <- loess.smooth(yvals, h0.est.func2, evaluation=500) ## Smooth out using loess
  
  ## Now using least squares to get estimated phi coefficients based on loess-smoothed h0
  h0.est.smooth <- h0.est.smooth.xy$y
  yvals.bp <- data.frame(predict(b, yvals))
  
  obj.fn.1 <- function(eta.vec, spline.pred, h0.truth){ 
    ##
    spline.vals <- sapply(1:nrow(spline.pred), function(i) sum(spline.pred[i,]*(eta.vec)))
    sq.diffs <- (log(h0.truth)-spline.vals)^2
    return(sum(sq.diffs))
  }
  if (method == "nlm"){
    phi.h0.est.smooth=nlm(f=obj.fn.1, p=rep(-1,num.BP.params), 
                          spline.pred=yvals.bp, 
                          h0.truth=h0.est.smooth)$estimate
    return(c(phi.h0.est.smooth, fitCox$coef))
  }
  if (method == "nls"){
    ## Alternatively -- this seems to give better results? 
    g <- function(bSpline.mat, phi.vec){
      phi.list <- as.list(phi.vec)
      mat.list <- as.list(data.frame(bSpline.mat))
      return(apply(do.call(rbind, Map('*', phi.list, mat.list)), 2, sum))
    }
    phi.h0.est.smooth.ls <- nls(h0.est.smooth ~ exp(g(yvals.bp, phi.vec)), 
                                start = list(phi.vec = rep(-1, ncol(yvals.bp))))
    phi.h0.est.smooth <- summary(phi.h0.est.smooth.ls)$parameters[,1]  ## These will be good for starting values!
    return(c(phi.h0.est.smooth, fitCox$coef))
  }
}


