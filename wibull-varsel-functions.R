#functions for BAR simulation study:
data.generator.12cov=function(turn,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv){
  set.seed(turn)
  
  if (grp.or.indiv=="indiv"){
    WB.simData=genWB.simData.noseed.allcontin.12cov (weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n =n,rho=rho, frailty = T,risks.status)
  }
  if (grp.or.indiv=="grp"){
    WB.simData=genWB.simData.noseed.allcontin.grp.12cov(G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n = n, rho=rho,frailty = T,risks.status)
  }
  
  p=(dim(WB.simData$Y)[2])-5
  #dimension when oracle:
  beta1.true.nonzero.length=sum(beta1.true!=0)
  beta2.true.nonzero.length=sum(beta2.true!=0)
  beta3.true.nonzero.length=sum(beta3.true!=0)
  p.oracle=sum(beta1.true.nonzero.length,beta2.true.nonzero.length,beta3.true.nonzero.length)
  
  #true beta values:
  beta.truth=(WB.simData$WBpara)[8:((3*p)+7)]
  #oracle true beta values:
  nonzero.param.truth.index=which((WB.simData$WBpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((WB.simData$WBpara)[8:((3*p)+7)])[nonzero.beta.truth.index]
  param.truth.oracle=(WB.simData$WBpara)[nonzero.param.truth.index]
  
  ## Data and model specifications
  Y <- WB.simData$Y                    # y1, delta1, y2, delta2, l = left-truncation time
  
  Y.oracle=(WB.simData$Y)[,-(5+(which(beta1.true==0)))]
  if (risks.status=="fullyshared"){
    Y.oracle.1=Y.oracle
    Y.oracle.2=Y.oracle
    Y.oracle.3=Y.oracle
  }
  if(risks.status=="overlapping"){
    Y.oracle.1=(WB.simData$Y)[,-(5+(which(beta1.true==0)))]
    Y.oracle.2=(WB.simData$Y)[,-(5+(which(beta2.true==0)))]
    Y.oracle.3=(WB.simData$Y)[,-(5+(which(beta3.true==0)))]
  }
  # 
  lin.pred <- WB.simData$lin.pred      # linear predictor formula
  
  lin.pred.oracle =WB.simData$lin.pred.oracle
  data <- WB.simData$data              # simulated data
  
  if (risks.status=="fullyshared"){
    data.oracle.1=data.oracle.2=data.oracle.3=WB.simData$data[,-(5+(which(beta1.true==0)))]
  }
  if (risks.status=="overlapping"){
    data.oracle.1=WB.simData$data[,-(5+(which(beta1.true==0)))]
    data.oracle.2=WB.simData$data[,-(5+(which(beta2.true==0)))]
    data.oracle.3=WB.simData$data[,-(5+(which(beta3.true==0)))]
    
  }
  # 
  pnum=length(beta1.true)
  trueWBparams <- WB.simData$WBpara    # true parameters used in data generation
  trueWBparams.oracle <- WB.simData$WBpara[nonzero.param.truth.index]
  frailty <- WB.simData$frailty        # allowing for a shared frailty
  para=WB.simData$WBpara
  para.oracle=trueWBparams.oracle 
  y1=WB.simData$Y[,1]
  y2=WB.simData$Y[,3]
  delta1=WB.simData$Y[,2]
  delta2=WB.simData$Y[,4]
  l=WB.simData$Y[,5]
  Ylength=dim(Y)[2]
  Ylength.oracle=dim(Y.oracle)[2]
  Y[,(6:Ylength)]->CovMat
  Y.oracle[,(6:Ylength.oracle)]->CovMat.oracle
  Xmat1=Xmat2=Xmat3=as.matrix(CovMat)
  if (risks.status=="fullyshared"){
    Xmat1.oracle=Xmat2.oracle=Xmat3.oracle=as.matrix(CovMat.oracle)
  }
  if(risks.status=="overlapping"){
    # Xmat1.oracle=as.matrix(CovMat[,nonzero.beta.truth.index[1:4]])
    # Xmat2.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[5:8])-10)])
    # Xmat3.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[9:12])-20)])
    Xmat1.oracle=as.matrix(CovMat[,nonzero.beta.truth.index[1:4]])
    Xmat2.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[5:8])-pnum)])
    Xmat3.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[9:12])-(2*pnum))])
    
  }
  n=dim(Xmat1)[1]
  E_ZZprime=(1/n)*(t(cbind(Xmat1,Xmat2,Xmat3))%*%cbind(Xmat1,Xmat2,Xmat3))
  
  p=ncol(Xmat3)
  p.oracle=ncol(Xmat3.oracle)
  # weibull.parameters=weibull.parameters
  # para.est=c(weibull.parameters,beta.truth)
  fit.wb=FreqID.LT(Y, lin.pred, data, model = "semi-Markov", frailty=frailty, startVals=NULL, method)
  
  if (risks.status=="fullyshared"){
    fit.wb.oracle=FreqID.LT(Y.oracle, lin.pred.oracle, data.oracle.1, model = "semi-Markov", frailty=frailty, startVals=NULL, method)
    
  }
  if (risks.status=="overlapping"){
    fit.wb.oracle=FreqID.LT.oracle.overlapping (Y.oracle.1, lin.pred.oracle,data.oracle.1,data.oracle.2,data.oracle.3, model = "semi-Markov", startVals=NULL, frailty=frailty,method)
    
  }
  
  
  
  # fit.wb.oracle=FreqID.LT(Y.oracle, lin.pred.oracle, data.oracle, model = "semi-Markov", frailty=frailty, startVals=NULL, method)
  #If cannot fit a model to find the unpenalized estimate, just jump over this iteration, and go to the next!
  if ('try-error' %in% class(fit.wb)){
    cat("I SEE A PROBLEM in running freqID.LT in data.generator function at seed:", seed," PLEASE CHECK THIS ITERATION \n")
    next
  }else{
    fit.wb=fit.wb
  }
  unpen.est=(fit.wb$estimate)[8:((3*p)+(8-1))]
  unpen.est.oracle=(fit.wb.oracle$estimate)[8:19]
  
  weibull.estimates=(fit.wb$estimate)[1:7]
  weibull.estimates.oracle=(fit.wb.oracle$estimate)[1:7]
  
  E_ZZprimei.oracle=list(0,n)
  
  for (i in 1:n){
    E_ZZprimei.oracle[[i]]=matrix(cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1)%*%t(matrix (cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1))
  }
  E_ZZprime.oracle=(1/n)*Reduce(`+`, E_ZZprimei.oracle)
  
  return(list(fit.wb=fit.wb,fit.wb.oracle=fit.wb.oracle,data=data,WB.simData=WB.simData,
              unpen.est.oracle=unpen.est.oracle, unpen.est=unpen.est,E_ZZprime=E_ZZprime,
              E_ZZprime.oracle=E_ZZprime.oracle,Y=Y,Y.oracle=Y.oracle,
              Y.oracle.1=Y.oracle.1,Y.oracle.2=Y.oracle.2,Y.oracle.3=Y.oracle.3,
              Xmat1.oracle=Xmat1.oracle,Xmat2.oracle=Xmat2.oracle,Xmat3.oracle=Xmat3.oracle,
              lin.pred=lin.pred,lin.pred.oracle=lin.pred.oracle,data=data,
              data.oracle.1=data.oracle.1,data.oracle.2=data.oracle.2,data.oracle.3=data.oracle.3,
              y1=y1,y2=y2,delta1=delta1,delta2=delta2,l=l,weibull.estimates=weibull.estimates,
              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,beta.truth.oracle=beta.truth.oracle,beta.truth=beta.truth))
}
#A function required i the modified shooting algorithm for Adaptiove Lasso Penalty Function
ss2 <- function(j,tmpb,Q,B)
{
  a <- sum(tmpb*Q[,j])-tmpb[j]*Q[j,j]
  s <- 2*(a-B[j])
  return(s)
}
#initi could be the unpenalized estimate: unpen.est
#Adaptive Lasso Penalty function solved using the modified shooting algorithm proposed in Zhang and Lu paper.
solveAdaLasso <- function(p,x,y,init,weight,lambda)
{
  Q = t(x)%*%x
  B = t(x)%*%y
  i=0
  status = 0
  
  lams =lambda*weight 
  oldbeta <- init
  tmpbeta <- oldbeta
  
  while (i<150 && status==0){ 
    for (j in 1:p){ 
      s <- ss2 (j,tmpbeta,Q,B)
      if (s > lams[j]) 
        tmpbeta[j]<-(lams[j]-s)/(2*Q[j,j])
      else if (s < (-lams[j]) ) 
        tmpbeta[j]<-(-lams[j]-s)/(2*Q[j,j])
      else
        tmpbeta[j]<- 0.0
    }
    dx<-max(abs(tmpbeta-oldbeta))
    oldbeta <- tmpbeta
    if (dx<=1e-06) 
      status <- 1 
    i <- i+1
  }
  return(tmpbeta)
}
#Lasso used Shooting algorithm to solve:
solveLasso <- function(yyy, XXX, lambda){  
  n=nrow(XXX)
  p=ncol(XXX)
  
  #Step 1: initialize beta, using reg. least square:
  XXXprime=t(XXX)
  first=XXXprime%*%XXX + 2*lambda 
  second=XXXprime%*%yyy
  # install.packages("pracma")
  library(pracma)
  beta <- mldivide(first,second)#vector of #23 elements
  
  #Step 2: for k=0,1,...,m, repeat:
  #convergence flag
  found <- 0
  # convergence tolerance
  TOL <- 1e-6
  while( found==0 ){
    #USing the current value of beta
    beta_old <- beta
    #optimize elements of beta by coordinate descent algorithm:
    for (i in 1:p){
      xxxi=XXX[,i]
      yyyi=yyy - XXX[,-i]%*%beta[-i,]
      xxxiprime=t(xxxi)
      deltai=xxxiprime%*%yyyi
      if (deltai< (-lambda)){
        firstt=deltai+lambda
        secondd=xxxiprime%*%xxxi
        beta[i]=firstt/secondd
      }
      else{
        if (deltai>lambda){
          firsttt=deltai-lambda
          seconddd=xxxiprime%*%xxxi
          beta[i]=firsttt/seconddd
        }
        else{
          beta[i]=0
        }
      }
      if (max(abs(beta-beta_old))<=TOL){
        found=1
      }
    }
  }
  #save outputs
  
  
  z.beta <- beta
  return(z.beta)
  
  
}
#BAR using the closed form solution to get solveD:
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
#The function used to solve bar with an iterative algorithm:
find.BARnew=function(lambdaopt,xiopt,para.est,y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3,frailty=TRUE) {
  nP.0 <- ifelse(frailty, 7, 6)
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  # beta= c(-0.03021803,-0.33064510,-0.10652843)  
  beta1=para.est[nP.0 + c(1:nP.1)]
  beta2=para.est[nP.0 + nP.1 + c(1:nP.2)]
  beta3=para.est[nP.0 + nP.1 + nP.2 + c(1:nP.3)]
  beta=c(beta1,beta2,beta3)
  flag=0
  TOL <- 1e-6
  count=0
  #fix other parameters other than beta:
  weibull.parameters=para.est[1:7]
  while (flag==0){
    betaold=beta
    weibull.parameters.old=weibull.parameters
    # print(betaold)
    #iteratively reweighted least square algorithm: 
    #first deriv. of lpglik function:
    para.est=c(weibull.parameters.old,betaold)
    G.ZhLunew=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    H.ZhLunew=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    # cat("para.est",para.est,"in iter:",count,"\n")
    X.ZhLunew=chol(H.ZhLunew)
    vecY.ZhLunew=forwardsolve(t(X.ZhLunew),H.ZhLunew%*%betaold-G.ZhLunew)
    beta=solveBAR(vecY.ZhLunew,X.ZhLunew,lambdaopt,xiopt)
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
  }
  # print(betaold) 
  # print(count)
  return(list(beta=beta,H=H.ZhLunew,G=G.ZhLunew,X=X.ZhLunew,y=vecY.ZhLunew))
} 
#finding ALASSO iteratively:
AdaLasso.finder.iter.solveAdaLasso=function(weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov", frailty=frailty, startVals=NULL,lam){
  para.est=c(weibull.estimates,unpen.est)
  G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(weibull.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(weibull.estimates,betaold)
    G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    # cat("para.est",para.est,"in iter:",count,"\n")
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    AdaLasso <-  solveAdaLasso((dim(X_irls)[1]),X_irls,Y_irls,unpen.est,1/abs(unpen.est),lam)
    # abar_cv <- cv.glmnet(x = X_irls, y = Y_irls,
    #                      type.measure = "mse",
    #                      nfold = 10,
    #                      alpha = 0,
    #                      penalty.factor = 1 / (abs(betaold))^2,
    #                      keep = TRUE,intercept=F)
    
    # best_abar_coef <- coef(abar_cv, s = lam)
    
    # if (count==1){
    #   print(abar)
    # }
    
    beta=AdaLasso
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  cat("last one:\n")
  # print(abar)
  MSE=t(as.matrix(beta-beta.truth))%*%(E_ZZprime)%*%as.matrix(beta-beta.truth)
  #number of true nonzero cov. :
  nonzero.truth=length(which(beta.truth!=0))
  # number of true zero cov. :
  zero.truth=length(which(beta.truth==0))
  nonzero.truth.index=which(beta.truth!=0)
  zero.truth.index=which(beta.truth==0)
  TP=length(which(abs(beta[nonzero.truth.index])!=0))
  FP=length(which(abs(beta[zero.truth.index])!=0))
  FN=length(which(abs(beta[nonzero.truth.index])==0))
  return(list(TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}
#Finding ALASSO iteratively in the group setting where covariates are clustered
AdaLasso.finder.iter.solveAdaLasso.grp=function(G1,G2,G3,G4,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov", frailty=frailty, startVals=NULL,lam){
  para.est=c(weibull.estimates,unpen.est)
  G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(weibull.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(weibull.estimates,betaold)
    G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    # cat("para.est",para.est,"in iter:",count,"\n")
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    AdaLasso <-  solveAdaLasso((dim(X_irls)[1]),X_irls,Y_irls,unpen.est,1/abs(unpen.est),lam)
    # abar_cv <- cv.glmnet(x = X_irls, y = Y_irls,
    #                      type.measure = "mse",
    #                      nfold = 10,
    #                      alpha = 0,
    #                      penalty.factor = 1 / (abs(betaold))^2,
    #                      keep = TRUE,intercept=F)
    
    # best_abar_coef <- coef(abar_cv, s = lam)
    
    # if (count==1){
    #   print(abar)
    # }
    
    beta=AdaLasso
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  cat("last one:\n")
  # print(abar)
  MSE=t(as.matrix(beta-beta.truth))%*%(E_ZZprime)%*%as.matrix(beta-beta.truth)
  #number of true nonzero cov. :
  nonzero.truth=length(which(beta.truth!=0))
  # number of true zero cov. :
  zero.truth=length(which(beta.truth==0))
  nonzero.truth.index=which(beta.truth!=0)
  zero.truth.index=which(beta.truth==0)
  TP=length(which(abs(beta[nonzero.truth.index])!=0))
  FP=length(which(abs(beta[zero.truth.index])!=0))
  # G1.extended=c(G1,10+G1,20+G1)
  # G2.extended=c(G2,10+G2,20+G2)
  # G3.extended=c(G3,10+G3,20+G3)
  # G4.extended=c(G4,10+G4,20+G4)
  
  pnum=(length(beta.truth))/3
  G1.extended=c(G1,pnum+G1,(2*pnum)+G1)
  G2.extended=c(G2,(pnum)+G2,(2*pnum)+G2)
  G3.extended=c(G3,pnum+G3,(2*pnum)+G3)
  G4.extended=c(G4,pnum+G4,(2*pnum)+G4)
  
  G1.perc=ifelse(sum(abs(beta[G1.extended])!=0)==length(G1.extended),1,0)
  G2.perc=ifelse(sum(abs(beta[G2.extended])!=0)==length(G2.extended),1,0)
  G3.perc=ifelse(sum(abs(beta[G3.extended])==0)==length(G3.extended),1,0)
  G4.perc=ifelse(sum(abs(beta[G4.extended])==0)==length(G4.extended),1,0)
  FN=length(which(abs(beta[nonzero.truth.index])==0))
  return(list(G1.perc=G1.perc,G2.perc=G2.perc,G3.perc=G3.perc,G4.perc=G4.perc,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}
#finding LASSO iteratively:
lasso.finder.iter.solvelasso=function(weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov", frailty=frailty, startVals=NULL,lam){
  para.est=c(weibull.estimates,unpen.est)
  G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(weibull.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(weibull.estimates,betaold)
    G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    # cat("para.est",para.est,"in iter:",count,"\n")
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    lasso <-  solveLasso(Y_irls,X_irls,lam)
    # abar_cv <- cv.glmnet(x = X_irls, y = Y_irls,
    #                      type.measure = "mse",
    #                      nfold = 10,
    #                      alpha = 0,
    #                      penalty.factor = 1 / (abs(betaold))^2,
    #                      keep = TRUE,intercept=F)
    
    # best_abar_coef <- coef(abar_cv, s = lam)
    
    # if (count==1){
    #   print(abar)
    # }
    
    beta=lasso
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  cat("last one:\n")
  # print(abar)
  MSE=t(as.matrix(beta-beta.truth))%*%(E_ZZprime)%*%as.matrix(beta-beta.truth)
  #number of true nonzero cov. :
  nonzero.truth=length(which(beta.truth!=0))
  # number of true zero cov. :
  zero.truth=length(which(beta.truth==0))
  nonzero.truth.index=which(beta.truth!=0)
  zero.truth.index=which(beta.truth==0)
  TP=length(which(beta[nonzero.truth.index]!=0))
  FP=length(which(beta[zero.truth.index]!=0))
  FN=length(which(beta[nonzero.truth.index]==0))
  return(list(TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}
#finding LASSO iteratively in group setting:
lasso.finder.iter.solvelasso.grp=function(G1,G2,G3,G4,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov", frailty=frailty, startVals=NULL,lam){
  para.est=c(weibull.estimates,unpen.est)
  G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(weibull.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(weibull.estimates,betaold)
    G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    # cat("para.est",para.est,"in iter:",count,"\n")
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    lasso <-  solveLasso(Y_irls,X_irls,lam)
    # abar_cv <- cv.glmnet(x = X_irls, y = Y_irls,
    #                      type.measure = "mse",
    #                      nfold = 10,
    #                      alpha = 0,
    #                      penalty.factor = 1 / (abs(betaold))^2,
    #                      keep = TRUE,intercept=F)
    
    # best_abar_coef <- coef(abar_cv, s = lam)
    
    # if (count==1){
    #   print(abar)
    # }
    
    beta=lasso
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  cat("last one:\n")
  # print(abar)
  MSE=t(as.matrix(beta-beta.truth))%*%(E_ZZprime)%*%as.matrix(beta-beta.truth)
  #number of true nonzero cov. :
  nonzero.truth=length(which(beta.truth!=0))
  # number of true zero cov. :
  zero.truth=length(which(beta.truth==0))
  nonzero.truth.index=which(beta.truth!=0)
  zero.truth.index=which(beta.truth==0)
  TP=length(which(beta[nonzero.truth.index]!=0))
  FP=length(which(beta[zero.truth.index]!=0))
  FN=length(which(beta[nonzero.truth.index]==0))
  
  pnum=(length(beta.truth))/3
  G1.extended=c(G1,pnum+G1,(2*pnum)+G1)
  G2.extended=c(G2,(pnum)+G2,(2*pnum)+G2)
  G3.extended=c(G3,pnum+G3,(2*pnum)+G3)
  G4.extended=c(G4,pnum+G4,(2*pnum)+G4)
  
  G1.perc=ifelse(sum(abs(beta[G1.extended])!=0)==length(G1.extended),1,0)
  G2.perc=ifelse(sum(abs(beta[G2.extended])!=0)==length(G2.extended),1,0)
  G3.perc=ifelse(sum(abs(beta[G3.extended])==0)==length(G3.extended),1,0)
  G4.perc=ifelse(sum(abs(beta[G4.extended])==0)==length(G4.extended),1,0)
  return(list(G1.perc=G1.perc,G2.perc=G2.perc,G3.perc=G3.perc,G4.perc=G4.perc,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}

BAR.finder.iter.solveBAR=function(tol,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov", frailty=frailty, startVals=NULL,lam,xi){
  para.est=c(weibull.estimates,unpen.est)
  G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(weibull.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(weibull.estimates,betaold)
    G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    # cat("para.est",para.est,"in iter:",count,"\n")
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    bar <-  solveBAR(Y_irls,X_irls,lam,xi)
    # abar_cv <- cv.glmnet(x = X_irls, y = Y_irls,
    #                      type.measure = "mse",
    #                      nfold = 10,
    #                      alpha = 0,
    #                      penalty.factor = 1 / (abs(betaold))^2,
    #                      keep = TRUE,intercept=F)
    
    # best_abar_coef <- coef(abar_cv, s = lam)
    
    # if (count==1){
    #   print(abar)
    # }
    
    beta=bar
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  cat("last one:\n")
  # print(abar)
  MSE=t(as.matrix(beta-beta.truth))%*%(E_ZZprime)%*%as.matrix(beta-beta.truth)
  #number of true nonzero cov. :
  nonzero.truth=length(which(beta.truth!=0))
  # number of true zero cov. :
  zero.truth=length(which(beta.truth==0))
  nonzero.truth.index=which(beta.truth!=0)
  zero.truth.index=which(beta.truth==0)
  TP=length(which(abs(beta[nonzero.truth.index])>tol))
  FP=length(which(abs(beta[zero.truth.index])>tol))
  FN=length(which(abs(beta[nonzero.truth.index])<=tol))
  return(list(TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}
BAR.finder.iter.solveBAR.grp=function(G1,G2,G3,G4,tol,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov", frailty=frailty, startVals=NULL,lam,xi){
  para.est=c(weibull.estimates,unpen.est)
  G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
  X_irls=chol(H)
  Y_irls=forwardsolve(t(X_irls),H%*%unpen.est-G)
  
  beta=unpen.est
  print(beta)
  # para.est=c(weibull.parameters,beta)
  para.est=c(weibull.estimates,beta)
  flag=0
  TOL <- 1e-6
  count=0
  while (flag==0 && count<=100){
    betaold=beta
    para.est=c(weibull.estimates,betaold)
    G=dlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    H=ddlogLike.weibull.new(para.est, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    # cat("para.est",para.est,"in iter:",count,"\n")
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    bar <-  solveBAR(Y_irls,X_irls,lam,xi)
    # abar_cv <- cv.glmnet(x = X_irls, y = Y_irls,
    #                      type.measure = "mse",
    #                      nfold = 10,
    #                      alpha = 0,
    #                      penalty.factor = 1 / (abs(betaold))^2,
    #                      keep = TRUE,intercept=F)
    
    # best_abar_coef <- coef(abar_cv, s = lam)
    
    # if (count==1){
    #   print(abar)
    # }
    
    beta=bar
    if (max(abs(beta-betaold))<=TOL){
      flag=1
    }
    count=count+1
    
    
  }
  cat("last one:\n")
  # print(abar)
  MSE=t(as.matrix(beta-beta.truth))%*%(E_ZZprime)%*%as.matrix(beta-beta.truth)
  #number of true nonzero cov. :
  nonzero.truth=length(which(beta.truth!=0))
  # number of true zero cov. :
  zero.truth=length(which(beta.truth==0))
  nonzero.truth.index=which(beta.truth!=0)
  zero.truth.index=which(beta.truth==0)
  TP=length(which(abs(beta[nonzero.truth.index])>tol))
  FP=length(which(abs(beta[zero.truth.index])>tol))
  FN=length(which(abs(beta[nonzero.truth.index])<=tol))
  
  # G1.extended=c(G1,10+G1,20+G1)
  # G2.extended=c(G2,10+G2,20+G2)
  # G3.extended=c(G3,10+G3,20+G3)
  # G4.extended=c(G4,10+G4,20+G4)
  pnum=(length(beta.truth))/3
  G1.extended=c(G1,pnum+G1,(2*pnum)+G1)
  G2.extended=c(G2,(pnum)+G2,(2*pnum)+G2)
  G3.extended=c(G3,pnum+G3,(2*pnum)+G3)
  G4.extended=c(G4,pnum+G4,(2*pnum)+G4)
  
  G1.perc=ifelse(sum(abs(beta[G1.extended])>tol)==length(G1.extended),1,0)
  G2.perc=ifelse(sum(abs(beta[G2.extended])>tol)==length(G2.extended),1,0)
  G3.perc=ifelse(sum(abs(beta[G3.extended])<=tol)==length(G3.extended),1,0)
  G4.perc=ifelse(sum(abs(beta[G4.extended])<=tol)==length(G4.extended),1,0)
  # G1.perc=ifelse(sum(abs(beta[G1])>tol)==length(G1),1,0)
  # G2.perc=ifelse(sum(abs(beta[G2])>tol)==length(G2),1,0)
  # G3.perc=ifelse(sum(abs(beta[G3])<=tol)==length(G3),1,0)
  # G4.perc=ifelse(sum(abs(beta[G4])<=tol)==length(G4),1,0)
  
  return(list(G1.perc=G1.perc,G2.perc=G2.perc,G3.perc=G3.perc,G4.perc=G4.perc,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}
#Finding lasso penalized estimate using generalized cross validetion criterion
GCV.finder.lasso=function(lambdavec,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL){
  for (lam in lambdavec){
    lasso.est=lasso.finder.iter.solvelasso(weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL,lam)
    TP.for.that.lam[which(lam==lambdavec)]=lasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=lasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=lasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=lasso.est$FN
    
    betavarsel.lasso[,(which(lambdavec==lam))]=lasso.est$beta
    betahat.and.weibull=c(weibull.estimates,lasso.est$beta)
    
    G.tilde.lasso=-ddlogLike.weibull.new(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
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
    numerator=logLike.weibull.SCR.SM.LT(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    denominator=n*(1-(p.lambda/n))^2
    GCV.lasso[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.lasso,"\n")
    optimal.gcv=which(GCV.lasso==min(GCV.lasso[!is.na(GCV.lasso)]))
    beta.GCV=betavarsel.lasso[,optimal.gcv]
  }
  GCV.lasso=GCV.lasso
  lam.final=lambdavec[optimal.gcv]
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  return(list(GCV.lasso=GCV.lasso,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.lasso,betavarsel.lasso=betavarsel.lasso))
}
#Finding lasso penalized estimate using generalized cross validetion criterion in group setting
GCV.finder.lasso.grp=function(G1,G2,G3,G4,lambdavec,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL){
  for (lam in lambdavec){
    lasso.est=lasso.finder.iter.solvelasso.grp(G1,G2,G3,G4,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL,lam)
    TP.for.that.lam[which(lam==lambdavec)]=lasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=lasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=lasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=lasso.est$FN
    G1.for.that.lam[which(lam==lambdavec)]=lasso.est$G1.perc
    G2.for.that.lam[which(lam==lambdavec)]=lasso.est$G2.perc
    G3.for.that.lam[which(lam==lambdavec)]=lasso.est$G3.perc
    G4.for.that.lam[which(lam==lambdavec)]=lasso.est$G4.perc
    
    betavarsel.lasso[,(which(lambdavec==lam))]=lasso.est$beta
    betahat.and.weibull=c(weibull.estimates,lasso.est$beta)
    
    G.tilde.lasso=-ddlogLike.weibull.new(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
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
    numerator=logLike.weibull.SCR.SM.LT(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    denominator=n*(1-(p.lambda/n))^2
    GCV.lasso[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.lasso,"\n")
    optimal.gcv=which(GCV.lasso==min(GCV.lasso[!is.na(GCV.lasso)]))
    beta.GCV=betavarsel.lasso[,optimal.gcv]
  }
  GCV.lasso=GCV.lasso
  lam.final=lambdavec[optimal.gcv]
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  pnum=(length(beta.truth))/3
  G1.extended=c(G1,pnum+G1,(2*pnum)+G1)
  G2.extended=c(G2,(pnum)+G2,(2*pnum)+G2)
  G3.extended=c(G3,pnum+G3,(2*pnum)+G3)
  G4.extended=c(G4,pnum+G4,(2*pnum)+G4)
  # G1.extended=c(G1,10+G1,20+G1)
  # G2.extended=c(G2,10+G2,20+G2)
  # G3.extended=c(G3,10+G3,20+G3)
  # G4.extended=c(G4,10+G4,20+G4)
  
  G1.perc.final=ifelse(sum(abs(beta.GCV[G1.extended])!=0)==length(G1.extended),1,0)
  G2.perc.final=ifelse(sum(abs(beta.GCV[G2.extended])!=0)==length(G2.extended),1,0)
  G3.perc.final=ifelse(sum(abs(beta.GCV[G3.extended])==0)==length(G3.extended),1,0)
  G4.perc.final=ifelse(sum(abs(beta.GCV[G4.extended])==0)==length(G4.extended),1,0)
  # G1.perc.final=G1.for.that.lam[optimal.gcv]
  # G2.perc.final=G2.for.that.lam[optimal.gcv]
  # G3.perc.final=G3.for.that.lam[optimal.gcv]
  # G4.perc.final=G4.for.that.lam[optimal.gcv]
  return(list(G1.perc.final=G1.perc.final,G2.perc.final=G2.perc.final,G3.perc.final=G3.perc.final,G4.perc.final=G4.perc.final,
              GCV.lasso=GCV.lasso,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.lasso,betavarsel.lasso=betavarsel.lasso))
}
#Finding Adaptive lasso penalized estimate using generalized cross validetion criterion
GCV.finder.AdaLasso=function(lambdavec,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL){
  for (lam in lambdavec){
    AdaLasso.est=AdaLasso.finder.iter.solveAdaLasso(weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL,lam)
    TP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FN
    
    betavarsel.AdaLasso[,(which(lambdavec==lam))]=AdaLasso.est$beta
    betahat.and.weibull=c(weibull.estimates,AdaLasso.est$beta)
    
    G.tilde.AdaLasso=-ddlogLike.weibull.new(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
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
    numerator=logLike.weibull.SCR.SM.LT(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    denominator=n*(1-(p.lambda/n))^2
    GCV.AdaLasso[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.AdaLasso,"\n")
    optimal.gcv=which(GCV.AdaLasso==min(GCV.AdaLasso[!is.na(GCV.AdaLasso)]))
    beta.GCV=betavarsel.AdaLasso[,optimal.gcv]
  }
  GCV.AdaLasso=GCV.AdaLasso
  lam.final=lambdavec[optimal.gcv]
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  return(list(GCV.AdaLasso=GCV.AdaLasso,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.AdaLasso,betavarsel.AdaLasso=betavarsel.AdaLasso))
}
#Finding Adaptive lasso penalized estimate using generalized cross validetion criterion in group setting
GCV.finder.AdaLasso.grp=function(G1,G2,G3,G4,lambdavec,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL){
  for (lam in lambdavec){
    AdaLasso.est=AdaLasso.finder.iter.solveAdaLasso.grp(G1,G2,G3,G4,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL,lam)
    TP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FN
    
    G1.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G1.perc
    G2.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G2.perc
    G3.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G3.perc
    G4.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G4.perc
    
    betavarsel.AdaLasso[,(which(lambdavec==lam))]=AdaLasso.est$beta
    betahat.and.weibull=c(weibull.estimates,AdaLasso.est$beta)
    
    G.tilde.AdaLasso=-ddlogLike.weibull.new(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
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
    numerator=logLike.weibull.SCR.SM.LT(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    denominator=n*(1-(p.lambda/n))^2
    GCV.AdaLasso[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.AdaLasso,"\n")
    optimal.gcv=which(GCV.AdaLasso==min(GCV.AdaLasso[!is.na(GCV.AdaLasso)]))
    beta.GCV=betavarsel.AdaLasso[,optimal.gcv]
  }
  GCV.AdaLasso=GCV.AdaLasso
  lam.final=lambdavec[optimal.gcv]
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  # G1.extended=c(G1,10+G1,20+G1)
  # G2.extended=c(G2,10+G2,20+G2)
  # G3.extended=c(G3,10+G3,20+G3)
  # G4.extended=c(G4,10+G4,20+G4)
  # 
  pnum=(length(beta.truth))/3
  G1.extended=c(G1,pnum+G1,(2*pnum)+G1)
  G2.extended=c(G2,(pnum)+G2,(2*pnum)+G2)
  G3.extended=c(G3,pnum+G3,(2*pnum)+G3)
  G4.extended=c(G4,pnum+G4,(2*pnum)+G4)
  
  G1.perc.final=ifelse(sum(abs(beta.GCV[G1.extended])!=0)==length(G1.extended),1,0)
  G2.perc.final=ifelse(sum(abs(beta.GCV[G2.extended])!=0)==length(G2.extended),1,0)
  G3.perc.final=ifelse(sum(abs(beta.GCV[G3.extended])==0)==length(G3.extended),1,0)
  G4.perc.final=ifelse(sum(abs(beta.GCV[G4.extended])==0)==length(G4.extended),1,0)
  # G1.perc.final=G1.for.that.lam[optimal.gcv]
  # G2.perc.final=G2.for.that.lam[optimal.gcv]
  # G3.perc.final=G3.for.that.lam[optimal.gcv]
  # G4.perc.final=G4.for.that.lam[optimal.gcv]
  return(list(G1.perc.final=G1.perc.final,G2.perc.final=G2.perc.final,G3.perc.final=G3.perc.final,G4.perc.final=G4.perc.final,
              GCV.AdaLasso=GCV.AdaLasso,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.AdaLasso,betavarsel.AdaLasso=betavarsel.AdaLasso))
}
#Finding BAR penalized estimate using generalized cross validetion criterion
GCV.finder.BAR=function(tol,lambdavec,xi,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL){
  for (lam in lambdavec){
    cat("this lam running:",lam,"\n")
    bar.est=BAR.finder.iter.solveBAR(tol,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL,lam,xi)
    
    TP.for.that.lam[which(lam==lambdavec)]=bar.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=bar.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=bar.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=bar.est$FN
    
    betavarsel.bar[,(which(lambdavec==lam))]=bar.est$beta
    betahat.and.weibull=c(weibull.estimates,bar.est$beta)
    
    G.tilde.bar=-ddlogLike.weibull.new(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
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
    numerator=logLike.weibull.SCR.SM.LT(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    denominator=n*(1-(p.lambda/n))^2
    GCV.bar[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.bar,"\n")
    optimal.gcv=which(GCV.bar==min(GCV.bar[!is.na(GCV.bar)]))
    beta.GCV=betavarsel.bar[,optimal.gcv]
  }
  GCV.bar=GCV.bar
  lam.final=lambdavec[optimal.gcv]
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  return(list(GCV.bar=GCV.bar,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.bar,betavarsel.bar=betavarsel.bar))
}
#Finding BAR penalized estimate using generalized cross validetion criterion in group setting
GCV.finder.BAR.grp=function(G1,G2,G3,G4,tol,lambdavec,xi,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL){
  for (lam in lambdavec){
    bar.est=BAR.finder.iter.solveBAR.grp(G1,G2,G3,G4,tol,weibull.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty=TRUE, startVals=NULL,lam,xi)
    TP.for.that.lam[which(lam==lambdavec)]=bar.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=bar.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=bar.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=bar.est$FN
    
    G1.for.that.lam[which(lam==lambdavec)]=bar.est$G1.perc
    G2.for.that.lam[which(lam==lambdavec)]=bar.est$G2.perc
    G3.for.that.lam[which(lam==lambdavec)]=bar.est$G3.perc
    G4.for.that.lam[which(lam==lambdavec)]=bar.est$G4.perc
    
    
    
    betavarsel.bar[,(which(lambdavec==lam))]=bar.est$beta
    betahat.and.weibull=c(weibull.estimates,bar.est$beta)
    
    G.tilde.bar=-ddlogLike.weibull.new(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
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
    numerator=logLike.weibull.SCR.SM.LT(betahat.and.weibull, y1, y2, delta1, delta2, l, Xmat1, Xmat2, Xmat3, frailty=TRUE)
    denominator=n*(1-(p.lambda/n))^2
    GCV.bar[which(lambdavec==lam)]=numerator/ denominator
    cat("this is GCV for lambda :",lam,"==>",GCV.bar,"\n")
    optimal.gcv=which(GCV.bar==min(GCV.bar[!is.na(GCV.bar)]))
    beta.GCV=betavarsel.bar[,optimal.gcv]
  }
  
  
  GCV.bar=GCV.bar
  lam.final=lambdavec[optimal.gcv]
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  # G1.extended=c(G1,10+G1,20+G1)
  # G2.extended=c(G2,10+G2,20+G2)
  # G3.extended=c(G3,10+G3,20+G3)
  # G4.extended=c(G4,10+G4,20+G4)
  
  pnum=(length(beta.truth))/3
  G1.extended=c(G1,pnum+G1,(2*pnum)+G1)
  G2.extended=c(G2,(pnum)+G2,(2*pnum)+G2)
  G3.extended=c(G3,pnum+G3,(2*pnum)+G3)
  G4.extended=c(G4,pnum+G4,(2*pnum)+G4)
  
  G1.perc.final=ifelse(sum(abs(beta.GCV[G1.extended])>tol)==length(G1.extended),1,0)
  G2.perc.final=ifelse(sum(abs(beta.GCV[G2.extended])>tol)==length(G2.extended),1,0)
  G3.perc.final=ifelse(sum(abs(beta.GCV[G3.extended])<=tol)==length(G3.extended),1,0)
  G4.perc.final=ifelse(sum(abs(beta.GCV[G4.extended])<=tol)==length(G4.extended),1,0)
  
  
  return(list(G1.perc.final=G1.perc.final,G2.perc.final=G2.perc.final,G3.perc.final=G3.perc.final,G4.perc.final=G4.perc.final,
              GCV.bar=GCV.bar,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.bar))
}

#Final result with replicates for some artificially generated data
repeat.get.last.result=function(seeds,tol,n,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,weibull.param.log,lambdavec.lasso,lambdavec.ada,lambdavec.bar,xi,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli){
  for (f in seeds){
    cat("this is risks.stat:",risks.status,"\n")
    cat("Now running seed",f,"\n")
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method)
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status)
    
    if(length(beta1.true)==10){
      dataa=try(data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==12){
      dataa=try(data.generator.12cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==14){
      dataa=try(data.generator.14cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==15){
      dataa=try(data.generator.15cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==16){
      dataa=try(data.generator.16cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==18){
      dataa=try(data.generator.15cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==19){
      dataa=try(data.generator.15cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==20){
      dataa=try(data.generator.20cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    
    unpen.est[which(f==seeds),]=dataa$unpen.est
    unpen.est.oracle[which(f==seeds),]=dataa$unpen.est.oracle
    # oracle.fit=FreqID.LT(dataa$Y.oracle, dataa$lin.pred.oracle, dataa$data.oracle, model = "semi-Markov", frailty=frailty, startVals=NULL)
    
    dimension[which(f==seeds)]=length(dataa$y1)
    
    censoring.observ.rates[which(f==seeds),1]=length(which(dataa$delta1==0&dataa$delta2==0))
    censoring.observ.rates[which(f==seeds),2]=length(which(dataa$delta1==0&dataa$delta2==1))
    censoring.observ.rates[which(f==seeds),3]=length(which(dataa$delta1==1&dataa$delta2==0))
    censoring.observ.rates[which(f==seeds),4]=length(which(dataa$delta1==1&dataa$delta2==1))
    
    # oracle.fit=try(FreqID.LT(dataa$Y.oracle.1, dataa$lin.pred.oracle, dataa$data.oracle.1, model = "semi-Markov", frailty=frailty, startVals=NULL,method),silent=TRUE)
    # if ('try-error' %in% class(oracle.fit)){
    #   cat("I see a problem in running FreqID.LT function, so, I am jumping to the next iteration \n")
    #   next
    # }else{
    #   oracle.fit=oracle.fit
    # }
    MSE.calculator.oracle[which(f==seeds)]=t(as.matrix((dataa$unpen.est.oracle)-dataa$beta.truth.oracle))%*%(dataa$E_ZZprime.oracle)%*%as.matrix((dataa$unpen.est.oracle)-dataa$beta.truth.oracle)
    # a=GCV.finder.BAR(tol,lambdavec,xi,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL)
    a=try(GCV.finder.BAR(tol,lambdavec.bar,xi,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(a)){
      cat("I see a problem in running GCV.finder function for BAR, so, I am jumping to the next iteration \n")
      next
    }else{
      a=a
    }
    blasso=try(GCV.finder.lasso(lambdavec.lasso,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(blasso)){
      cat("I see a problem in running GCV.finder function for lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      blasso=blasso
    }
    cada=try(GCV.finder.AdaLasso(lambdavec.ada,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    
    # a=try(GCV.finder.lasso(lambdavec,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(cada)){
      cat("I see a problem in running GCV.finder function for adaptive lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      cada=cada
    }
    # GCV.selected[(which(seed==seeds))]=a$lam.final
    GCV.selected.bar[,(which(f==seeds))]=a$GCV.bar
    selected.beta.based.on.GCV.bar=a$beta.GCV
    TP.bar[(which(f==seeds))]=a$TP.final
    FP.bar[(which(f==seeds))]=a$FP.final
    MSE.bar[(which(f==seeds))]=a$MSE.final
    beta.selected.bar[,(which(f==seeds))]=a$beta.GCV
    FP.bar[(which(f==seeds))]=a$FP.final
    FN.bar[(which(f==seeds))]=a$FN.final
    
    GCV.selected.lasso[,(which(f==seeds))]=blasso$GCV.lasso
    selected.beta.based.on.GCV.lasso=blasso$beta.GCV
    TP.lasso[(which(f==seeds))]=blasso$TP.final
    beta.selected.lasso[,(which(f==seeds))]=blasso$beta.GCV
    FP.lasso[(which(f==seeds))]=blasso$FP.final
    MSE.lasso[(which(f==seeds))]=blasso$MSE.final
    FN.lasso[(which(f==seeds))]=blasso$FN.final
    
    GCV.selected.ada[,(which(f==seeds))]=cada$GCV.AdaLasso
    selected.beta.based.on.GCV.ada=cada$beta.GCV
    TP.ada[(which(f==seeds))]=cada$TP.final
    beta.selected.ada[,(which(f==seeds))]=cada$beta.GCV
    FP.ada[(which(f==seeds))]=cada$FP.final
    MSE.ada[(which(f==seeds))]=cada$MSE.final
    FN.ada[(which(f==seeds))]=cada$FN.final
    
    
  }
  
  censoring.rate.final=apply(censoring.observ.rates,2,mean)
  dimension.final.mean=mean(dimension[!is.na(dimension)])
  
  n.worked.lasso=length(seeds)-sum(is.na(TP.lasso))
  n.worked.ada=length(seeds)-sum(is.na(TP.ada))
  n.worked.bar=length(seeds)-sum(is.na(TP.bar))
  n.worked.oracle=length(seeds)-sum(is.na(MSE.calculator.oracle))
  
  SD.lasso=sd(MSE.lasso[!is.na(MSE.lasso)])*sqrt((n.worked.lasso-1)/(n.worked.lasso))
  SD.ada=sd(MSE.ada[!is.na(MSE.ada)])*sqrt((n.worked.ada-1)/(n.worked.ada))
  SD.bar=sd(MSE.bar[!is.na(MSE.bar)])*sqrt((n.worked.bar-1)/(n.worked.bar))
  SD.oracle=sd(MSE.calculator.oracle[!is.na(MSE.calculator.oracle)])*sqrt((n.worked.oracle-1)/(n.worked.oracle))
  
  MC.lasso=FP.lasso+FN.lasso
  MC.ada=FP.ada+FN.ada
  MC.bar=FP.bar+FN.bar
  
  MC.final.lasso=mean(MC.lasso[!is.na(MC.lasso)])
  MC.final.ada=mean(MC.ada[!is.na(MC.ada)])
  MC.final.bar=mean(MC.bar[!is.na(MC.bar)])
  
  return(list(censoring.rate.final=censoring.rate.final,censoring.observ.rates=censoring.observ.rates,dimension=dimension,dimension.final.mean=dimension.final.mean,
              unpen.est=unpen.est,unpen.est.oracle=unpen.est.oracle
              ,SD.bar=SD.bar,beta.selected.bar=beta.selected.bar,GCV.selected.bar=GCV.selected.bar,TP.bar=TP.bar,FP.bar=FP.bar,FN.bar=FN.bar,MC.final.bar=MC.final.bar,MSE.bar=MSE.bar,MSE.calculator.oracle=MSE.calculator.oracle
              ,SD.lasso=SD.lasso,beta.selected.lasso=beta.selected.lasso,GCV.selected.lasso=GCV.selected.lasso,TP.lasso=TP.lasso,FP.lasso=FP.lasso,FN.lasso=FN.lasso,MC.lasso=MC.lasso,MC.final.lasso=MC.final.lasso,MSE.lasso=MSE.lasso
              ,SD.ada=SD.ada,beta.selected.ada=beta.selected.ada,GCV.selected.ada=GCV.selected.ada,TP.ada=TP.ada,FP.ada=FP.ada,FN.ada=FN.ada,MC.ada=MC.ada,MC.final.ada=MC.final.ada,MSE.ada=MSE.ada,SD.oracle=SD.oracle))
}
##Final result with replicates for some artificially grouped simulated data
repeat.get.last.result.grp=function(seeds,tol,n,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,weibull.param.log,lambdavec.lasso,lambdavec.ada,lambdavec.bar,xi,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli){
  for (f in seeds){
    cat("this is risks.stat:",risks.status,"\n")
    cat("Now running seed",f,"\n")
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv)
    
    
    if(length(beta1.true)==12){
      dataa=try(data.generator.12cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    
    if(length(beta1.true)==15){
      dataa=try(data.generator.15cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if(length(beta1.true)==16){
      dataa=try(data.generator.16cov(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    
    # dataa=try(data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv),silent=TRUE)
    # if ('try-error' %in% class(dataa)){
    #   cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
    #   next
    # }else{
    #   dataa=dataa
    # }
    # oracle.fit=FreqID.LT(dataa$Y.oracle, dataa$lin.pred.oracle, dataa$data.oracle, model = "semi-Markov", frailty=frailty, startVals=NULL)
    
    dimension[which(f==seeds)]=length(dataa$y1)
    
    censoring.observ.rates[which(f==seeds),1]=length(which(dataa$delta1==0&dataa$delta2==0))
    censoring.observ.rates[which(f==seeds),2]=length(which(dataa$delta1==0&dataa$delta2==1))
    censoring.observ.rates[which(f==seeds),3]=length(which(dataa$delta1==1&dataa$delta2==0))
    censoring.observ.rates[which(f==seeds),4]=length(which(dataa$delta1==1&dataa$delta2==1))
    
    # oracle.fit=try(FreqID.LT(dataa$Y.oracle, dataa$lin.pred.oracle, dataa$data.oracle, model = "semi-Markov", frailty=frailty, startVals=NULL,method),silent=TRUE)
    # if ('try-error' %in% class(oracle.fit)){
    #   cat("I see a problem in running FreqID.LT function, so, I am jumping to the next iteration \n")
    #   next
    # }else{
    #   oracle.fit=oracle.fit
    # }
    MSE.calculator.oracle[which(f==seeds)]=t(as.matrix((dataa$unpen.est.oracle)-dataa$beta.truth.oracle))%*%(dataa$E_ZZprime.oracle)%*%as.matrix((dataa$unpen.est.oracle)-dataa$beta.truth.oracle)
    # a=GCV.finder.BAR.grp(G1,G2,G3,G4,tol,lambdavec.bar,xi,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL)
    
    a=try(GCV.finder.BAR.grp(G1,G2,G3,G4,tol,lambdavec.bar,xi,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(a)){
      cat("I see a problem in running GCV.finder function, so, I am jumping to the next iteration \n")
      next
    }else{
      a=a
    }
    
    # blasso=GCV.finder.lasso.grp(G1,G2,G3,G4,lambdavec.lasso,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL)
    
    blasso=try(GCV.finder.lasso.grp(G1,G2,G3,G4,lambdavec.lasso,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(blasso)){
      cat("I see a problem in running GCV.finder function for lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      blasso=blasso
    }
    
    # cada=GCV.finder.AdaLasso.grp(G1,G2,G3,G4,lambdavec.ada,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL)
    
    cada=try(GCV.finder.AdaLasso.grp(G1,G2,G3,G4,lambdavec.ada,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    
    # a=try(GCV.finder.lasso(lambdavec,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(cada)){
      cat("I see a problem in running GCV.finder function for adaptive lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      cada=cada
    }
    # GCV.selected[(which(seed==seeds))]=a$lam.final
    GCV.selected.bar[,(which(f==seeds))]=a$GCV.bar
    selected.beta.based.on.GCV.bar=a$beta.GCV
    TP.bar[(which(f==seeds))]=a$TP.final
    FP.bar[(which(f==seeds))]=a$FP.final
    MSE.bar[(which(f==seeds))]=a$MSE.final
    beta.selected.bar[,(which(f==seeds))]=a$beta.GCV
    FN.bar[(which(f==seeds))]=a$FN.final
    G1.percentage.bar[(which(f==seeds))]=a$G1.perc.final
    G2.percentage.bar[(which(f==seeds))]=a$G2.perc.final
    G3.percentage.bar[(which(f==seeds))]=a$G3.perc.final
    G4.percentage.bar[(which(f==seeds))]=a$G4.perc.final
    
    GCV.selected.lasso[,(which(f==seeds))]=blasso$GCV.lasso
    selected.beta.based.on.GCV.lasso=blasso$beta.GCV
    TP.lasso[(which(f==seeds))]=blasso$TP.final
    beta.selected.lasso[,(which(f==seeds))]=blasso$beta.GCV
    FP.lasso[(which(f==seeds))]=blasso$FP.final
    MSE.lasso[(which(f==seeds))]=blasso$MSE.final
    FN.lasso[(which(f==seeds))]=blasso$FN.final
    G1.percentage.lasso[(which(f==seeds))]=blasso$G1.perc.final
    G2.percentage.lasso[(which(f==seeds))]=blasso$G2.perc.final
    G3.percentage.lasso[(which(f==seeds))]=blasso$G3.perc.final
    G4.percentage.lasso[(which(f==seeds))]=blasso$G4.perc.final
    
    GCV.selected.ada[,(which(f==seeds))]=cada$GCV.AdaLasso
    selected.beta.based.on.GCV.ada=cada$beta.GCV
    TP.ada[(which(f==seeds))]=cada$TP.final
    beta.selected.ada[,(which(f==seeds))]=cada$beta.GCV
    FP.ada[(which(f==seeds))]=cada$FP.final
    MSE.ada[(which(f==seeds))]=cada$MSE.final
    FN.ada[(which(f==seeds))]=cada$FN.final
    G1.percentage.ada[(which(f==seeds))]=cada$G1.perc.final
    G2.percentage.ada[(which(f==seeds))]=cada$G2.perc.final
    G3.percentage.ada[(which(f==seeds))]=cada$G3.perc.final
    G4.percentage.ada[(which(f==seeds))]=cada$G4.perc.final
    
    
  }
  
  censoring.rate.final=apply(censoring.observ.rates,2,mean)
  
  n.worked.lasso=length(seeds)-sum(is.na(TP.lasso))
  n.worked.ada=length(seeds)-sum(is.na(TP.ada))
  n.worked.bar=length(seeds)-sum(is.na(TP.bar))
  n.worked.oracle=length(seeds)-sum(is.na(MSE.calculator.oracle))
  
  SD.lasso=sd(MSE.lasso[!is.na(MSE.lasso)])*sqrt((n.worked.lasso-1)/(n.worked.lasso))
  SD.ada=sd(MSE.ada[!is.na(MSE.ada)])*sqrt((n.worked.ada-1)/(n.worked.ada))
  SD.bar=sd(MSE.bar[!is.na(MSE.bar)])*sqrt((n.worked.bar-1)/(n.worked.bar))
  SD.oracle=sd(MSE.calculator.oracle[!is.na(MSE.calculator.oracle)])*sqrt((n.worked.oracle-1)/(n.worked.oracle))
  
  MC.lasso=FP.lasso+FN.lasso
  MC.ada=FP.ada+FN.ada
  MC.bar=FP.bar+FN.bar
  
  MC.final.lasso=mean(MC.lasso[!is.na(MC.lasso)])
  MC.final.ada=mean(MC.ada[!is.na(MC.ada)])
  MC.final.bar=mean(MC.bar[!is.na(MC.bar)])
  
  G1.percentage.f.bar=sum(G1.percentage.bar[!is.na(G1.percentage.bar)])/length(G1.percentage.bar[!is.na(G1.percentage.bar)])
  G2.percentage.f.bar=sum(G2.percentage.bar[!is.na(G2.percentage.bar)])/length(G2.percentage.bar[!is.na(G2.percentage.bar)])
  G3.percentage.f.bar=sum(G3.percentage.bar[!is.na(G3.percentage.bar)])/length(G3.percentage.bar[!is.na(G3.percentage.bar)])
  G4.percentage.f.bar=sum(G4.percentage.bar[!is.na(G4.percentage.bar)])/length(G4.percentage.bar[!is.na(G4.percentage.bar)])
  # G.total.bar=(0.2*G1.percentage.f.bar)+(0.2*G2.percentage.f.bar)+(0.3*G3.percentage.f.bar)+(0.3*G4.percentage.f.bar)
  G.total.bar=((length(G1)/length(beta1.true))*G1.percentage.f.bar)+((length(G2)/length(beta1.true))*G2.percentage.f.bar)+((length(G3)/length(beta1.true))*G3.percentage.f.bar)+((length(G4)/length(beta1.true))*G4.percentage.f.bar)
  
  
  G1.percentage.f.lasso=sum(G1.percentage.lasso[!is.na(G1.percentage.lasso)])/length(G1.percentage.lasso[!is.na(G1.percentage.lasso)])
  G2.percentage.f.lasso=sum(G2.percentage.lasso[!is.na(G2.percentage.lasso)])/length(G2.percentage.lasso[!is.na(G2.percentage.lasso)])
  G3.percentage.f.lasso=sum(G3.percentage.lasso[!is.na(G3.percentage.lasso)])/length(G3.percentage.lasso[!is.na(G3.percentage.lasso)])
  G4.percentage.f.lasso=sum(G4.percentage.lasso[!is.na(G4.percentage.lasso)])/length(G4.percentage.lasso[!is.na(G4.percentage.lasso)])
  # G.total.lasso=(0.2*G1.percentage.f.lasso)+(0.2*G2.percentage.f.lasso)+(0.3*G3.percentage.f.lasso)+(0.3*G4.percentage.f.lasso)
  G.total.lasso=((length(G1)/length(beta1.true))*G1.percentage.f.lasso)+((length(G2)/length(beta1.true))*G2.percentage.f.lasso)+((length(G3)/length(beta1.true))*G3.percentage.f.lasso)+((length(G4)/length(beta1.true))*G4.percentage.f.lasso)
  
  
  G1.percentage.f.ada=sum(G1.percentage.ada[!is.na(G1.percentage.ada)])/length(G1.percentage.ada[!is.na(G1.percentage.ada)])
  G2.percentage.f.ada=sum(G2.percentage.ada[!is.na(G2.percentage.ada)])/length(G2.percentage.ada[!is.na(G2.percentage.ada)])
  G3.percentage.f.ada=sum(G3.percentage.ada[!is.na(G3.percentage.ada)])/length(G3.percentage.ada[!is.na(G3.percentage.ada)])
  G4.percentage.f.ada=sum(G4.percentage.ada[!is.na(G4.percentage.ada)])/length(G4.percentage.ada[!is.na(G4.percentage.ada)])
  # G.total.ada=(0.2*G1.percentage.f.ada)+(0.2*G2.percentage.f.ada)+(0.3*G3.percentage.f.ada)+(0.3*G4.percentage.f.ada)
  # 
  G.total.ada=((length(G1)/length(beta1.true))*G1.percentage.f.ada)+((length(G2)/length(beta1.true))*G2.percentage.f.ada)+((length(G3)/length(beta1.true))*G3.percentage.f.ada)+((length(G4)/length(beta1.true))*G4.percentage.f.ada)
  
  return(list(G.total.bar=G.total.bar,G1.percentage.f.bar=G1.percentage.f.bar,G2.percentage.f.bar=G2.percentage.f.bar,G3.percentage.f.bar=G3.percentage.f.bar,G4.percentage.f.bar=G4.percentage.f.bar,
              G.total.lass=G.total.lasso,G1.percentage.f.lasso=G1.percentage.f.lasso,G2.percentage.f.lasso=G2.percentage.f.lasso,G3.percentage.f.lasso=G3.percentage.f.lasso,G4.percentage.f.lasso=G4.percentage.f.lasso,
              G.total.ada=G.total.ada,G1.percentage.f.ada=G1.percentage.f.ada,G2.percentage.f.ada=G2.percentage.f.ada,G3.percentage.f.ada=G3.percentage.f.ada,G4.percentage.f.ada=G4.percentage.f.ada,
              censoring.rate.final=censoring.rate.final,censoring.observ.rates=censoring.observ.rates,dimension=dimension
              ,SD.bar=SD.bar,beta.selected.bar=beta.selected.bar,GCV.selected.bar=GCV.selected.bar,TP.bar=TP.bar,FP.bar=FP.bar,FN.bar=FN.bar,MC.final.bar=MC.final.bar,MSE.bar=MSE.bar,MSE.calculator.oracle=MSE.calculator.oracle
              ,SD.lasso=SD.lasso,beta.selected.lasso=beta.selected.lasso,GCV.selected.lasso=GCV.selected.lasso,TP.lasso=TP.lasso,FP.lasso=FP.lasso,FN.lasso=FN.lasso,MC.lasso=MC.lasso,MC.final.lasso=MC.final.lasso,MSE.lasso=MSE.lasso
              ,SD.ada=SD.ada,beta.selected.ada=beta.selected.ada,GCV.selected.ada=GCV.selected.ada,TP.ada=TP.ada,FP.ada=FP.ada,FN.ada=FN.ada,MC.ada=MC.ada,MC.final.ada=MC.final.ada,MSE.ada=MSE.ada,SD.oracle=SD.oracle))
}
