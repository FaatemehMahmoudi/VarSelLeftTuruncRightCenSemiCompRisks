#########Application of project 2##############
#rm(list=ls())
#for Bernstein Polynomials:
library(dplyr)
library(survival)
require(stats)
library(splines2)
library(pracma)  ## for numerical differentiation
library(MASS)
# setwd("/Users/fatemehmahmoudi/Desktop/Codes-bin")

source("sim-functions-bp.R")
source("varsel-functions-bp.R")
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

FreqID.LT.real.data <- function(Y, lin.pred, data, model = "semi-Markov", startVals, frailty=TRUE, method="nlm")
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
  fit.survreg.1 <- survreg(as.formula(paste("Surv(y1, delta1) ", as.character(lin.pred[[1]])[1], as.character(lin.pred[[1]])[2])), dist="weibull", data=data)
  fit.survreg.2 <- survreg(as.formula(paste("Surv(y2, delta2) ", as.character(lin.pred[[2]])[1], as.character(lin.pred[[2]])[2])), dist="weibull", data=data)
  # data.delta1_1 = data[delta1==1,]
  # data.delta1_1$y2.m.y1 = y2[delta1==1] - y1[delta1==1]
  data.delta1_1 = data[delta1==1,]
  data.delta1_1$y2.m.y1 = y2[delta1==1] - y1[delta1==1]
  data.delta1_1=data.delta1_1[-c(which(data.delta1_1$y2.m.y1==0)),]
  fit.survreg.3 <- survreg(as.formula(paste("Surv(y2.m.y1, delta2) ", as.character(lin.pred[[3]])[1], as.character(lin.pred[[3]])[2])), dist="weibull", data=data.delta1_1)
  alpha1      <- 1 / fit.survreg.1$scale
  alpha2      <- 1 / fit.survreg.2$scale
  alpha3     	<- 1 / fit.survreg.3$scale
  
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
  
  value <- list(estimate=est, H=H, logLike=logLike, myLabels=myLabels, frailty=frailty, nP=nP, 
                code = code)
  
  class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
  
  return(value)
  ##
  invisible()
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

solveBAR <- function(Y, X, lambda,xi){  
  
  p=ncol(X)
  #Step 1: initialize beta, using ridge reg:
  Im <- diag(1, p, p)
  beta <- solve(t(X)%*% X + xi*Im)%*%t(X)%*%Y
  #Step 2: for k=0,1,...,m, repeat:
  #convergence flag
  found <- 0
  # convergence tolerance
  TOL <- 1e-6
  count=0
  d <- 0.001 #to prevent computation overflow:
  while( found==0 ){
    #USing the current value of beta
    beta_old <- beta
    D <- diag( as.vector(1/(beta_old^2 + d^2)), p, p)
    beta <- solve(t(X)%*%X + lambda*D)%*%t(X)%*%Y
    count=count+1
    if (count>100){
      break;
    }
    if (max(abs(beta-beta_old))<=TOL){
      found=1
    }
  }
  #save outputs
  z.beta <- beta
  return(z.beta)
}

AdaLasso.finder.iter.solveAdaLasso.bp=function(nuisance.estimates,unpen.est,Y,y1,y2,
                                               delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred, data, 
                                               model = "semi-Markov",frailty=frailty, startVals,lam,
                                               b.1,b.2,b.3.y2my1){
  para.est=c(nuisance.estimates,unpen.est)
  G=dlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                b.1,  
                                                b.2,  
                                                b.3.y2my1)
  H=ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                 b.1,  
                                                 b.2,  
                                                 b.3.y2my1)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  para.est=c(nuisance.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(nuisance.estimates,betaold)
    G=dlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                  b.1,
                                                  b.2,
                                                  b.3.y2my1)
    H=ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                   b.1,
                                                   b.2,
                                                   b.3.y2my1)
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    AdaLasso <-  solveAdaLasso(3*p,X_irls,Y_irls,unpen.est,1/abs(unpen.est),lam)

    beta=AdaLasso
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1


  }
  return(list(beta=beta,H=H,G=G,X=X_irls,y=Y_irls,count=count))
  
}

lasso.finder.iter.solvelasso.bp=function(nuisance.estimates,unpen.est,Y,y1,y2,delta1,
                                         delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
                                         frailty=frailty, startVals,lam,
                                         b.1,b.2,b.3.y2my1){
  para.est=c(nuisance.estimates,unpen.est)
  G=dlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2,l, Xmat1, Xmat2, Xmat3, frailty,
                                                b.1,  
                                                b.2,  
                                                b.3.y2my1)
  H=ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                 b.1,  
                                                 b.2,  
                                                 b.3.y2my1)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(nuisance.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(nuisance.estimates,betaold)
    G=dlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                  b.1,  
                                                  b.2,  
                                                  b.3.y2my1)
    H=ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                   b.1,  
                                                   b.2,  
                                                   b.3.y2my1)
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    lasso <-  solveLasso(Y_irls,X_irls,lam)
   
    beta=lasso
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  cat("last one:\n")
  # print(abar)
  
  return(list(beta=beta,H=H,G=G,X=X_irls,y=Y_irls,count=count))
  
}

bar.finder.iter.solvebar.bp=function(tol,nuisance.estimates,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
                                     frailty=frailty,lam,xi,b.1,b.2,b.3.y2my1){
  
  para.est=c(nuisance.estimates,unpen.est)
  G=dlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE,
                                                b.1,  
                                                b.2,  
                                                b.3.y2my1)
  H=ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE,
                                                 b.1,  
                                                 b.2,  
                                                 b.3.y2my1)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(nuisance.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(nuisance.estimates,betaold)
    G=dlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE,
                                                  b.1,  
                                                  b.2,  
                                                  b.3.y2my1)
    H=ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE,
                                                   b.1,  
                                                   b.2,  
                                                   b.3.y2my1)
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    bar <-  solveBAR(Y_irls,X_irls,lam,xi)
    
    beta=bar
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  
  return(list(beta=beta,H=H,G=G,X=X_irls,y=Y_irls,count=count))
  
}
GCV.finder.AdaLasso.bp=function(lambdavec,nuisance.estimates,unpen.est,Y,y1,y2,delta1,
                                delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty,
                                startVals,
                                b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    cat("running lam=",lam,"\n")
    AdaLasso.est=AdaLasso.finder.iter.solveAdaLasso.bp(nuisance.estimates,unpen.est,Y,y1,y2,
                                                       delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, 
                                                       model = "semi-Markov", frailty=frailty, startVals,lam,
                                                       b.1,b.2,b.3.y2my1)
    
    
    betavarsel.AdaLasso[,(which(lambdavec==lam))]=AdaLasso.est$beta
    betahat.and.nuisance=c(nuisance.estimates,AdaLasso.est$beta)
    
    G.tilde.AdaLasso=-ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(betahat.and.nuisance, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                                   b.1,  
                                                                   b.2,  
                                                                   b.3.y2my1)
    s.beta=c(rep(NA,length(AdaLasso.est$beta)))
    for (i in 1:length(AdaLasso.est$beta)){
      if (AdaLasso.est$beta[i]==0){
        s.beta[i]=0.000001
      }else{
        s.beta[i]=abs(AdaLasso.est$beta[i])
      }
    }
    n=dim(Xmat1)[1]
    A=diag(1/s.beta)
    first=solve(G.tilde.AdaLasso+(lam*A))
    
    p.lambda=sum(diag(first%*%G.tilde.AdaLasso))
    numerator=logLike.SCR.SM.LT.bSpline.bp.dropPrevCases(betahat.and.nuisance, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                         b.1,  
                                                         b.2,  
                                                         b.3.y2my1)
    denominator=n*(1-(p.lambda/n))^2
    GCV.AdaLasso[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.AdaLasso,"\n")
    optimal.gcv=which(GCV.AdaLasso==min(GCV.AdaLasso[!is.na(GCV.AdaLasso)]))
    beta.GCV=betavarsel.AdaLasso[,optimal.gcv]
  }
  GCV.AdaLasso=GCV.AdaLasso
  lam.final=lambdavec[optimal.gcv]
  
  
  return(list(GCV.AdaLasso=GCV.AdaLasso,lam.final=lam.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.AdaLasso,betavarsel.AdaLasso=betavarsel.AdaLasso))
}
GCV.finder.lasso.bp=function(lambdavec,nuisance.estimates,unpen.est,Y,y1,y2,delta1,
                             delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov",
                             frailty, startVals,
                             b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    lasso.est=lasso.finder.iter.solvelasso.bp(nuisance.estimates,unpen.est,Y,y1,y2,delta1,
                                              delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
                                              frailty, startVals,lam,
                                              b.1,b.2,b.3.y2my1)
    
    betavarsel.lasso[,(which(lambdavec==lam))]=lasso.est$beta
    betahat.and.nuisance=c(nuisance.estimates,lasso.est$beta)
    
    G.tilde.lasso=-ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(betahat.and.nuisance, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                                b.1,  
                                                                b.2,  
                                                                b.3.y2my1)
    s.beta=c(rep(NA,length(lasso.est$beta)))
    for (i in 1:length(lasso.est$beta)){
      if (lasso.est$beta[i]==0){
        s.beta[i]=0.000001
      }else{
        s.beta[i]=abs(lasso.est$beta[i])
      }
    }
    n=dim(Xmat1)[1]
    A=diag(1/s.beta)
    first=solve(G.tilde.lasso+(lam*A))
    
    p.lambda=sum(diag(first%*%G.tilde.lasso))
    
    numerator=logLike.SCR.SM.LT.bSpline.bp.dropPrevCases(betahat.and.nuisance, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                         b.1,  
                                                         b.2,  
                                                         b.3.y2my1)
    denominator=n*(1-(p.lambda/n))^2
    GCV.lasso[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.lasso,"\n")
    optimal.gcv=which(GCV.lasso==min(GCV.lasso[!is.na(GCV.lasso)]))
    beta.GCV=betavarsel.lasso[,optimal.gcv]
  }
  GCV.lasso=GCV.lasso
  lam.final=lambdavec[optimal.gcv]
  
  return(list(GCV.lasso=GCV.lasso,lam.final=lam.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.lasso,betavarsel.lasso=betavarsel.lasso))
}

GCV.finder.BAR.bp=function(tol,lambdavec,xi,nuisance.estimates,unpen.est,Y,y1,y2,
                           delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov",
                           frailty, startVals,
                           b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    cat("this lam running:",lam,"\n")
    bar.est=bar.finder.iter.solvebar.bp(tol,nuisance.estimates,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
                                        frailty=frailty,lam,xi,b.1,b.2,b.3.y2my1)
    
    betavarsel.bar[,(which(lambdavec==lam))]=bar.est$beta
    betahat.and.nuisance=c(nuisance.estimates,bar.est$beta)
    
    G.tilde.bar=-ddlogLike.SCR.SM.LT.bSpline.bp.dropPrevCases(betahat.and.nuisance, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                              b.1,  
                                                              b.2,  
                                                              b.3.y2my1)      
    s.beta=c(rep(NA,length(bar.est$beta)))
    for (i in 1:length(bar.est$beta)){
      if (bar.est$beta[i]==0){
        s.beta[i]=0.000001
      }else{
        s.beta[i]=abs(bar.est$beta[i])
      }
    }
    n=dim(Xmat1)[1]
    A=diag(1/s.beta)
    first=solve(G.tilde.bar+(lam*A))
    
    p.lambda=sum(diag(first%*%G.tilde.bar))
    numerator=logLike.SCR.SM.LT.bSpline.bp.dropPrevCases(betahat.and.nuisance, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty,
                                                         b.1,  
                                                         b.2,  
                                                         b.3.y2my1)    
    denominator=n*(1-(p.lambda/n))^2
    GCV.bar[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.bar,"\n")
    optimal.gcv=which(GCV.bar==min(GCV.bar[!is.na(GCV.bar)]))
    beta.GCV=betavarsel.bar[,optimal.gcv]
  }
  GCV.bar=GCV.bar
  lam.final=lambdavec[optimal.gcv]
  return(list(GCV.bar=GCV.bar,lam.final=lam.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.bar,betavarsel.bar=betavarsel.bar))
}

code_start_time<-Sys.time()
data<-colon
data<-data[!is.na(data[,9]),]
data<-data[!is.na(data[,11]),]
num_type=2
data_length<-length(data$id)
n<-data_length/2
maxiter=200
eps=1e-4
p=12
K=2
num_tuning=20
replicate=1
r=1
beta_hat_lasso_array<-array(0,dim=c(p,K,replicate))
beta_hat_lasso_adap_array<-array(0,dim=c(p,K,replicate))
beta_hat_origin_array<-array(0,dim=c(p,K,replicate))
beta_hat_array<-array(0,dim=c(p,K,replicate))
beta_hat_adap_array<-array(0,dim=c(p,K,replicate))
adapweight<-array(0,dim=c(p,K,replicate))


Data_analysis<-function(start,t,delta,x,adapweight_bi_level,tuning,maxiter,eps=1e-4,trace=FALSE){
  
  meanx<-matrix(0,p,K)
  normx<-matrix(0,p,K)
  covariate<-array(0,dim=c(n,p,K))
  x.gamma<-array(0,dim=c(n,p,K))
  beta.scaled<-matrix(0,p,K)
  for(k in 1:K){
    covariate[,,k] = scale(x[,,k])
    meanx[,k] = attributes(scale(x[,,k]))$'scaled:center'
    normx[,k] = attributes(scale(x[,,k]))$'scaled:scale'
  }
  
  
  ## initial value ##
  gamma = rep(1,p)
  gamma.old = gamma
  theta<-matrix(0,p,K)
  theta.old = matrix(0,p,K)
  if(missing(adapweight_bi_level)){
    adapweight_bi_level=matrix(1,p,K)
  }
  iter_record<-0
  iter = 0
  dif = 1
  
  while(iter<maxiter && dif>eps) {
    iter = iter + 1
    ## update theta ##
    for(k in 1:K){
      x.theta=t(t(covariate[,,k])*gamma)
      x.theta = scale(x.theta, FALSE, adapweight_bi_level[,k]) ## incorporate adaptive weights
      fit.theta = penalized(Surv(start,t[,k],delta[,k]), penalized=x.theta, unpenalized = ~0, model='cox', lambda1=tuning, standardize=FALSE)
      theta[,k] = as.vector(attributes(fit.theta)$penalized)/adapweight_bi_level[,k]
      x.gamma[,,k]= t(t(covariate[,,k])*theta[,k])
    }
    
    ## Update gamma Using Cycle Coordinate Descent##
    maxiter_gamma=maxiter
    eps_gamma=0.5
    iter_gamma=0
    dif_gamma=1
    gamma_origin=gamma
    gamma_update=gamma
    u_gamma<-rep(0,p)
    A_gamma<-rep(0,p)
    while(iter_gamma<maxiter_gamma && dif_gamma>eps_gamma){
      iter_gamma=iter_gamma+1
      for(j in 1:p){
        if(theta[j,1]==0&theta[j,2]==0){
          gamma_update[j]=0
          next
        }else{
          for(i in 1:n){
            for(k in 1:K){
              S_1=0
              S=0
              S_2=0
              for(m in 1:n){
                betax_l=gamma_origin%*%x.gamma[m,,k]
                S_1=ifelse(t[m,k]>=t[i,k],1,0)*exp(betax_l)*x.gamma[m,j,k]+S_1
                S=ifelse(t[m,k]>=t[i,k],1,0)*exp(betax_l)+S
                S_2=ifelse(t[m,k]>=t[i,k],1,0)*exp(betax_l)*((x.gamma[m,j,k])^2)+S_2
              }
              u_gamma[j]=u_gamma[j]+delta[i,k]*(x.gamma[i,j,k]-S_1/S)
              A_gamma[j]=A_gamma[j]+delta[i,k]*(-S_2/S+(S_1/S)^2)
            }  
          }
          u_gamma[j]=u_gamma[j]
          A_gamma[j]=A_gamma[j]
          if(gamma_origin[j]>=0){
            gamma_update[j]=max(0,gamma_origin[j]-(u_gamma[j]-tuning)/A_gamma[j])
          }else{
            gamma_update[j]=min(0,gamma_origin[j]-(u_gamma[j]+tuning)/A_gamma[j])
          }
        }
      }
      dif_gamma=max(abs(gamma_update-gamma_origin))
      gamma_origin=gamma_update
      cat('iter_gamma:', iter_gamma, '\n')
      cat('dif_gamma:', dif_gamma, '\n')
    } 
    gamma= as.vector(abs(gamma_origin)) 
    
    dif = max(max(abs(theta-theta.old)), max(abs(gamma-gamma.old)))    
    theta.old = theta
    gamma.old = gamma
    
    cat('iter:', iter, '\n')
    cat('dif:', dif, '\n')
  }
  
  beta.scaled= theta*gamma
  beta_hat =beta.scaled/normx
  
  iter_record<-iter  
  return(beta_hat)
}


t=matrix(0,data_length/num_type,num_type)
delta=matrix(0,data_length/num_type,num_type)
count1_t=1
count2_t=1
count1_delta=1
count2_delta=1
etype<-data$etype
for(i in 1:data_length){
  if(etype[i]==1){
    t[count1_t,1]=data$time[i]
    delta[count1_delta,1]=data$status[i]
    count1_t=count1_t+1 
    count1_delta=count1_delta+1  
  } else {
    
    t[count2_t,2]=data$time[i]
    delta[count2_delta,2]=data$status[i]
    count2_t=count2_t+1
    count2_delta=count2_delta+1
  }
}

######Generating Covariates#####
x<-array(0,dim=c(data_length/2,12,num_type))
rx<-data$rx
count1_lev=1
count2_lev=1
count1_lev5fu=1
count2_lev5fu=1
for(i in 1:data_length){
  if(etype[i]==1){
    if(rx[i]=="Lev"){
      x[count1_lev,1,1]=1
    } else {
      x[count1_lev,1,1]=0
    }
    if(rx[i]=="Lev+5FU"){
      x[count1_lev,2,1]=1
    } else {
      x[count1_lev,2,1]=0
    }
    count1_lev=count1_lev+1
    count1_lev5fu=count1_lev5fu+1
  } else {
    if(rx[i]=="Lev"){
      x[count1_lev,1,2]=1
    } else {
      x[count1_lev,1,2]=0
    }
    if(rx[i]=="Lev+5FU"){
      x[count1_lev,2,2]=1
    } else {
      x[count1_lev,2,2]=0
    }
    count2_lev=count2_lev+1
    count2_lev5fu=count2_lev5fu+1
  }
} 

count1_sex=1
count2_sex=1
count1_age=1
count2_age=1
count1_obs=1
count2_obs=1
count1_per=1
count2_per=1
count1_adh=1
count2_adh=1
count1_nods=1
count2_nods=1
count1_dif=1
count2_dif=1
count1_ext=1
count2_ext=1
count1_sur=1
count2_sur=1
count1_nod4=1
count2_nod4=1
etype<-data$etype
for(i in 1:data_length){
  if(etype[i]==1){
    x[count1_sex,3,1]=data$sex[i]
    x[count1_age,4,1]=data$age[i]
    x[count1_obs,5,1]=data$obstruct[i]
    x[count1_per,6,1]=data$perfor[i]
    x[count1_adh,7,1]=data$adhere[i]
    x[count1_nods,8,1]=data$nodes[i]
    x[count1_dif,9,1]=data$differ[i]
    x[count1_ext,10,1]=data$extent[i]
    x[count1_sur,11,1]=data$surg[i]
    x[count1_nod4,12,1]=data$node4[i]
    count1_sex=count1_sex+1 
    count1_age=count1_age+1
    count1_obs=count1_obs+1 
    count1_per=count1_per+1
    count1_adh=count1_adh+1 
    count1_nods=count1_nods+1
    count1_dif=count1_dif+1 
    count1_ext=count1_ext+1
    count1_sur=count1_sur+1 
    count1_nod4=count1_nod4+1
  } else {
    x[count2_sex,3,2]=data$sex[i]
    x[count2_age,4,2]=data$age[i]
    x[count2_obs,5,2]=data$obstruct[i]
    x[count2_per,6,2]=data$perfor[i]
    x[count2_adh,7,2]=data$adhere[i]
    x[count2_nods,8,2]=data$nodes[i]
    x[count2_dif,9,2]=data$differ[i]
    x[count2_ext,10,2]=data$extent[i]
    x[count2_sur,11,2]=data$surg[i]
    x[count2_nod4,12,2]=data$node4[i]
    count2_sex=count2_sex+1 
    count2_age=count2_age+1
    count2_obs=count2_obs+1 
    count2_per=count2_per+1
    count2_adh=count2_adh+1 
    count2_nods=count2_nods+1
    count2_dif=count2_dif+1 
    count2_ext=count2_ext+1
    count2_sur=count2_sur+1 
    count2_nod4=count2_nod4+1
  }
}


meanx<-matrix(0,p,K)
normx<-matrix(0,p,K)
covariate<-array(0,dim=c(n,p,K))
for(k in 1:K){
  covariate[,,k] = scale(x[,,k])
  meanx[,k] = attributes(scale(x[,,k]))$'scaled:center'
  normx[,k] = attributes(scale(x[,,k]))$'scaled:scale'
}

Zcov=x[,,1]
Zcov=as.matrix(Zcov)


Zcovscale=covariate[,,1]
Zcovscale=as.matrix(Zcovscale)

start <- rep(0,n)
#first transition:
stop=t[,1]
time1=stop-start
status1=delta[,1]

#second transition:
stop=t[,2]
time2=stop-start
status2=delta[,2]


firstportionofdata=cbind(time1,status1,time2,status2)

dim(firstportionofdata)[1]==dim(Zcov)[1]
dim(firstportionofdata)[1]==dim(Zcovscale)[1]

#creating an artificial vector of left truncations all being zero:
L.crc=c(rep(0,dim(firstportionofdata)[1]))
dat.crc.reg=cbind(firstportionofdata,L.crc,Zcov)
dat.crc.scaled=cbind(firstportionofdata,L.crc,Zcovscale)

colnames(dat.crc.scaled)[1:17]=c("y1","delta1","y2","delta2","L","first","second","third","fourth","fifth","sixth","seventh","eighth","ninth","tenth","eleventh","twelfth")
colnames(dat.crc.reg)[1:17]=c("y1","delta1","y2","delta2","L","first","second","third","fourth","fifth","sixth","seventh","eighth","ninth","tenth","eleventh","twelfth")

lin.pred.crc= list(as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth), as.formula(~first+second+third+fourth+fifth+sixth+seventh+eighth+ninth+tenth+eleventh+twelfth))

dat.crc.reg=as.data.frame(dat.crc.reg)
dat.crc.scaled=as.data.frame(dat.crc.scaled)

#some data are wrongly entered so that y2-y1 for those for whom delta1 and delta2 are both 1 are 
#the same. It should not be the case in semi-comp data. So, I dlete those cases which are 4 cases:

dat.crc.reg$del=dat.crc.reg$y2-dat.crc.reg$y1
indextodelete=which(dat.crc.reg$del==0&dat.crc.reg$delta1==1&dat.crc.reg$delta2==1)

dat.crc.reg=dat.crc.reg[-c(indextodelete),]
dat.crc.reg=dat.crc.reg[,-18]

dat.crc.scaled=dat.crc.scaled[-c(indextodelete),]
dat.crc.scaled=dat.crc.scaled[,-18]

dim(dat.crc.reg)

#unpenalized estimate with BP BH:
m1=9
m2=9
m3=9
y1=dat.crc.reg$y1
y2=dat.crc.reg$y2
delta1=dat.crc.reg$delta1
delta2=dat.crc.reg$delta2
lin.pred=lin.pred.crc
data=dat.crc.reg
#in the bernstein polynomial method, we use the function "approxfun" in r to approximate the 
#hazard function based on a set of points on x axis (l for example), and a SET OF POINTS ON Y AXIS (bernstein basis values )
#to build the terms in likelihood function while in Weibull method we do not need this procedure. 
#this "approxfun" function does not work if all the values on x axis are the same. So, here 
#I am artificially produce some very small values almost near zero to handle this issue:

Y=data
data$L=sample(seq(0,0.001,0.0001),dim(Y)[1],replace=TRUE)
Y=data
head(data)
bdy.knots.b.1 = c(0, max(y1)) 
bdy.knots.b.2 = c(0, max(y2))
bdy.knots.b.3.y2my1 = c(0, max(y2-y1))

#first trans:
b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
b.1.bp=predict(b.1.event.bp,y1)
start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength=500, method="nls")
phi.1.bp.start=start1.bp[1:(m1+1)]
beta1.bp.start=start1.bp[(m1+1+1):(dim(Zcov)[2]+(m1+1))]
#second trans:
b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
b.2.bp=predict(b.2.event.bp,y2)
start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength=500, method="nls")
phi.2.bp.start=start2.bp[1:(m2+1)]
beta2.bp.start=start2.bp[(m2+1+1):(dim(Zcov)[2]+(m2+1))]
#third trans:
b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength=500, method="nls")
phi.3.bp.start=start3.bp[1:(m3+1)]
beta3.bp.start=start3.bp[(m3+1+1):(dim(Zcov)[2]+(m3+1))]

startVals.bp=c(phi.1.bp.start,phi.2.bp.start,phi.3.bp.start,log(0.25),beta1.bp.start,beta2.bp.start,beta3.bp.start)

#without scaling the covariates:
fit.bh.crc.reg=FreqID.LT.bSpline.bp(Y,lin.pred, data,startVals=startVals.bp, frailty=TRUE, 
                     b.1=b.1.bp, b.2=b.2.bp, b.3.y2my1=b.3.y2my1.bp, 
                     bdy.knots.b.1, 
                     bdy.knots.b.2, 
                     bdy.knots.b.3.y2my1,
                     method="optim")

#in the bernstein polynomial method, we use the function "approxfun" in r to approximate the 
#hazard function based on a set of points on x axis (l for example), and a SET OF POINTS ON Y AXIS (bernstein basis values )
#to build the terms in likelihood function while in Weibull method we do not need this procedure. 
#this "approxfun" function does not work if all the values on x axis are the same. So, here 
#I am artificially produce some very small values almost near zero to handle this issue:

Y=data
data$L=sample(seq(0,0.001,0.0001),dim(Y)[1],replace=TRUE)
Y=data
head(data)
bdy.knots.b.1 = c(0, max(y1)) 
bdy.knots.b.2 = c(0, max(y2))
bdy.knots.b.3.y2my1 = c(0, max(y2-y1))

#first trans:
b.1.event.bp=bernsteinPoly(x=y1[delta1==1],degree=m1,intercept = TRUE,Boundary.knots = bdy.knots.b.1)
b.1.bp=predict(b.1.event.bp,y1)
start1.bp=get.BP.startVals(m1,y1, delta1, lin.pred[[1]], data, b.1.bp,yvalslength=500, method="nls")
phi.1.bp.start=start1.bp[1:(m1+1)]
beta1.bp.start=start1.bp[(m1+1+1):(dim(Zcov)[2]+(m1+1))]
#second trans:
b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
b.2.bp=predict(b.2.event.bp,y2)
start2.bp=get.BP.startVals(m2,y2, delta2, lin.pred[[2]], data, b.2.bp,yvalslength=500, method="nls")
phi.2.bp.start=start2.bp[1:(m2+1)]
beta2.bp.start=start2.bp[(m2+1+1):(dim(Zcov)[2]+(m2+1))]
#third trans:
b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
start3.bp=get.BP.startVals(m3,(y2-y1)[delta1==1], delta2[delta1==1], lin.pred[[3]], data[delta1==1,],b.3.y2my1.bp,yvalslength=500, method="nls")
phi.3.bp.start=start3.bp[1:(m3+1)]
beta3.bp.start=start3.bp[(m3+1+1):(dim(Zcov)[2]+(m3+1))]

startVals.bp=c(phi.1.bp.start,phi.2.bp.start,phi.3.bp.start,log(0.25),beta1.bp.start,beta2.bp.start,beta3.bp.start)

fit.bh.reg=FreqID.LT.bSpline.bp(Y,lin.pred, data,startVals=startVals.bp, frailty=TRUE, 
                                    b.1=b.1.bp, b.2=b.2.bp, b.3.y2my1=b.3.y2my1.bp, 
                                    bdy.knots.b.1, 
                                    bdy.knots.b.2, 
                                    bdy.knots.b.3.y2my1,
                                    method="optim")


#censoring rates:
censoring.observ.rates=matrix(NA,1,4)

censoring.observ.rates[1]=length(which(data$delta1==0&data$delta2==0))
censoring.observ.rates[2]=length(which(data$delta1==0&data$delta2==1))
censoring.observ.rates[3]=length(which(data$delta1==1&data$delta2==0))
censoring.observ.rates[4]=length(which(data$delta1==1&data$delta2==1))

censoring.observ.rates=as.data.frame(censoring.observ.rates)
colnames(censoring.observ.rates)=c("00","01","10","11")
censoring.observ.rates


###################################################
######      VARIABLE SELECTION     ################
###################################################


#data with scaled version of covariates:

data=dat.crc.scaled
Y=data
data$L=sample(seq(0,0.001,0.0001),dim(Y)[1],replace=TRUE)
Y=data
l=data$L
head(data)
# lambdavec.bar=seq(1.86,2.48,0.05)
# lambdavec.lasso=seq(5.5,6.5,0.05)
# lambdavec.ada=seq(5.5,6.5,0.05)
lambdavec.bar=seq(3,3.1,0.05)
lambdavec.lasso=seq(5.5,6.5,0.05)
lambdavec.ada=seq(5.5,6.5,0.05)
tol=1e-04
xi=10
nuisance.estimates.s=fit.bh.crc.scaled$estimate[1:(m1+1+m2+1+m3+1+1)]
unpen.est.s=fit.bh.crc.scaled$estimate[(m1+1+m2+1+m3+1+1+1):length(fit.bh.crc.scaled$estimate)]
Xmat1.s=Xmat2.s=Xmat3.s=as.matrix(data[,c(6:dim(data)[2])])
p=12
GCV.bar=vector()
GCV.selected.bar=matrix(0,length(lambdavec.bar),1)
betavarsel.bar=matrix(0,3*p,length(lambdavec.bar))
beta.selected.bar=matrix(0,3*p,1)
GCV.lasso=vector()
GCV.selected.lasso=matrix(0,length(lambdavec.lasso),1)
betavarsel.lasso=matrix(0,3*p,length(lambdavec.lasso))
beta.selected.lasso=matrix(0,3*p,1)
GCV.AdaLasso=vector()
GCV.selected.ada=matrix(0,length(lambdavec.ada),1)
betavarsel.AdaLasso=matrix(0,3*p,length(lambdavec.ada))
beta.selected.ada=matrix(0,3*p,1)

#BAR:
a.bp.bar=GCV.finder.BAR.bp(tol,lambdavec.bar,xi,nuisance.estimates.s,unpen.est,Y,
                              data$y1,data$y2,data$delta1,data$delta2,data$L,Xmat1.s,Xmat2.s,
                              Xmat3.s,lin.pred,data,model = "semi-Markov",
                              frailty=TRUE,startVals.bp,b.1.bp ,b.2.bp,b.3.y2my1.bp)
a.scale$betavarsel.bar->crc.bar.scale.res
crc.bar.scale.res=as.data.frame(crc.bar.scale.res)
rownames(crc.bar.scale.res)=c("Lev1","Lev+FU1","Sex1","Age1","Obstruct1","Perfor1","Adhere1","Nodes1","Differ1","Extent1","Surg1","Node41","Lev2","Lev+FU2","Sex2","Age2","Obstruct2","Perfor2","Adhere2","Nodes2","Differ2","Extent2","Surg2","Node42","Lev3","Lev+FU3","Sex3","Age3","Obstruct3","Perfor3","Adhere3","Nodes3","Differ3","Extent3","Surg3","Node413")
colnames(crc.bar.scale.res)=lambdavec.bar
Answer.scaled.bar.GCV=a.scale$beta.GCV
Answer.scaled.bar.GCV=as.data.frame(Answer.scaled.bar.GCV)
rownames(Answer.scaled.bar.GCV)=c("Lev1","Lev+FU1","Sex1","Age1","Obstruct1","Perfor1","Adhere1","Nodes1","Differ1","Extent1","Surg1","Node41","Lev2","Lev+FU2","Sex2","Age2","Obstruct2","Perfor2","Adhere2","Nodes2","Differ2","Extent2","Surg2","Node42","Lev3","Lev+FU3","Sex3","Age3","Obstruct3","Perfor3","Adhere3","Nodes3","Differ3","Extent3","Surg3","Node413")
indextobezeroinbar=which(abs(Answer.scaled.bar.GCV)<=tol)
Answer.scaled.bar.GCV[indextobezeroinbar,]=0
#LASSO:
b.scaled=GCV.finder.lasso.bp(lambdavec.lasso,nuisance.estimates,unpen.est,Y,data$y1,data$y2,data$delta1,
                             data$delta2,data$L,Xmat1,Xmat2,Xmat3,lin.pred,data,model = "semi-Markov",
                             frailty=TRUE, startVals.bp,
                             b.1.bp,b.2.bp,b.3.y2my1.bp)
b.scaled$betavarsel.lasso->crc.lasso.scale.res
crc.lasso.scale.res=as.data.frame(crc.lasso.scale.res)
rownames(crc.lasso.scale.res)=c("Lev1","Lev+FU1","Sex1","Age1","Obstruct1","Perfor1","Adhere1","Nodes1","Differ1","Extent1","Surg1","Node41","Lev2","Lev+FU2","Sex2","Age2","Obstruct2","Perfor2","Adhere2","Nodes2","Differ2","Extent2","Surg2","Node42","Lev3","Lev+FU3","Sex3","Age3","Obstruct3","Perfor3","Adhere3","Nodes3","Differ3","Extent3","Surg3","Node413")
colnames(crc.lasso.scale.res)=lambdavec.lasso
Answer.scaled.lasso.GCV=b.scaled$beta.GCV
Answer.scaled.lasso.GCV=as.data.frame(Answer.scaled.lasso.GCV)
rownames(Answer.scaled.lasso.GCV)=c("Lev1","Lev+FU1","Sex1","Age1","Obstruct1","Perfor1","Adhere1","Nodes1","Differ1","Extent1","Surg1","Node41","Lev2","Lev+FU2","Sex2","Age2","Obstruct2","Perfor2","Adhere2","Nodes2","Differ2","Extent2","Surg2","Node42","Lev3","Lev+FU3","Sex3","Age3","Obstruct3","Perfor3","Adhere3","Nodes3","Differ3","Extent3","Surg3","Node413")

#ADAPTIVE LASSO:
c.scaled=GCV.finder.AdaLasso.bp(lambdavec.ada,nuisance.estimates,unpen.est,Y,data$y1,data$y2,data$delta1,
                                data$delta2,data$L,Xmat1,Xmat2,Xmat3,lin.pred,data,model = "semi-Markov",
                                frailty=TRUE,startVals.bp,
                                b.1.bp,b.2.bp,b.3.y2my1.bp)
c.scaled$betavarsel.AdaLasso->crc.adalasso.scale.res
crc.adalasso.scale.res=as.data.frame(crc.adalasso.scale.res)
rownames(crc.adalasso.scale.res)=c("Lev1","Lev+FU1","Sex1","Age1","Obstruct1","Perfor1","Adhere1","Nodes1","Differ1","Extent1","Surg1","Node41","Lev2","Lev+FU2","Sex2","Age2","Obstruct2","Perfor2","Adhere2","Nodes2","Differ2","Extent2","Surg2","Node42","Lev3","Lev+FU3","Sex3","Age3","Obstruct3","Perfor3","Adhere3","Nodes3","Differ3","Extent3","Surg3","Node413")
colnames(crc.adalasso.scale.res)=lambdavec.lasso
Answer.scaled.adalasso.GCV=c.scaled$beta.GCV
Answer.scaled.adalasso.GCV=as.data.frame(Answer.scaled.adalasso.GCV)
rownames(Answer.scaled.adalasso.GCV)=c("Lev1","Lev+FU1","Sex1","Age1","Obstruct1","Perfor1","Adhere1","Nodes1","Differ1","Extent1","Surg1","Node41","Lev2","Lev+FU2","Sex2","Age2","Obstruct2","Perfor2","Adhere2","Nodes2","Differ2","Extent2","Surg2","Node42","Lev3","Lev+FU3","Sex3","Age3","Obstruct3","Perfor3","Adhere3","Nodes3","Differ3","Extent3","Surg3","Node413")


final.result=cbind(unpen.est,Answer.scaled.bar.GCV,Answer.scaled.adalasso.GCV,Answer.scaled.lasso.GCV)
print(final.result)

setwd("/Users/fatemehmahmoudi/Desktop/Codes-bin/real-data-crc-jun6")
getwd()
filename=c("CRC-BAR-LASS-ALASS-UNP")
sink(filename)
print(final.result)
sink()






