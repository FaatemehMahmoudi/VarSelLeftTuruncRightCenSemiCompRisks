data.generator.and.unpenalized.estimator=function(turn,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution){
  set.seed(turn)
  # # beta1.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta2.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta3.true=c(-1,-1,-1,-1,0,0,0,0,0,0)
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  
  
  if(grp.or.indiv=="indiv"){
    BP.simData=genWB.simData.10cov.BP.BS(mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution)
  }
  if (grp.or.indiv=="grp"){
    BP.simData=genWB.simData.10cov.BP.BS.grp(G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,n ,m,yvalslength, rho,frailty ,risks.status)
  }
  
  p=(dim(BP.simData$Y)[2])-5
  #dimension when oracle:
  beta1.true.nonzero.length=sum(beta1.true!=0)
  beta2.true.nonzero.length=sum(beta2.true!=0)
  beta3.true.nonzero.length=sum(beta3.true!=0)
  p.oracle=sum(beta1.true.nonzero.length,beta2.true.nonzero.length,beta3.true.nonzero.length)
  
  #true beta values:
  nuisance.p= sum(BP.simData$m)+3+1
  beta.truth=(BP.simData$BSpara)[(nuisance.p+1):((3*p)+nuisance.p)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$BSpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$BPpara)[(nuisance.p+1):((3*p)+nuisance.p)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$BPpara)[nonzero.param.truth.index]
  ## Data and model specifications
  Y <- BP.simData$Y                    # y1, delta1, y2, delta2, l = left-truncation time
  
  Y.oracle=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
  if (risks.status=="fullyshared"){
    Y.oracle.1=Y.oracle
    Y.oracle.2=Y.oracle
    Y.oracle.3=Y.oracle
  }
  if(risks.status=="overlapping"){
    Y.oracle.1=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
    Y.oracle.2=(BP.simData$Y)[,-(5+(which(beta2.true==0)))]
    Y.oracle.3=(BP.simData$Y)[,-(5+(which(beta3.true==0)))]
  }
  
  lin.pred <- BP.simData$lin.pred      # linear predictor formula
  lin.pred.oracle =BP.simData$lin.pred.oracle
  
  data <- BP.simData$data              # simulated data
  
  if (risks.status=="fullyshared"){
    data.oracle.1=data.oracle.2=data.oracle.3=BP.simData$data[,-(5+(which(beta1.true==0)))]
  }
  if (risks.status=="overlapping"){
    data.oracle.1=BP.simData$data[,-(5+(which(beta1.true==0)))]
    data.oracle.2=BP.simData$data[,-(5+(which(beta2.true==0)))]
    data.oracle.3=BP.simData$data[,-(5+(which(beta3.true==0)))]
    
  }
  
  trueBPparams <- BP.simData$BPpara    # true parameters used in data generation
  trueBPparams.oracle <- BP.simData$BPpara[nonzero.param.truth.index]
  frailty <- BP.simData$frailty        # allowing for a shared frailty
  para=BP.simData$BPpara
  para.oracle=trueBPparams.oracle 
  y1=BP.simData$Y[,1]
  y2=BP.simData$Y[,3]
  delta1=BP.simData$Y[,2]
  delta2=BP.simData$Y[,4]
  l=BP.simData$Y[,5]
  Ylength=dim(Y)[2]
  Ylength.oracle=dim(Y.oracle)[2]
  Y[,(6:Ylength)]->CovMat
  Y.oracle[,(6:Ylength.oracle)]->CovMat.oracle
  
  
  Xmat1=Xmat2=Xmat3=as.matrix(CovMat)
  if (risks.status=="fullyshared"){
    Xmat1.oracle=Xmat2.oracle=Xmat3.oracle=as.matrix(CovMat.oracle)
  }
  p=ncol(Xmat1)
  if(risks.status=="overlapping"){
    Xmat1.oracle=as.matrix(CovMat[,nonzero.beta.truth.index[1:4]])
    Xmat2.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[5:8])-p)])
    Xmat3.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[9:12])-(2*p))])
    
  }
  
  
  E_ZZprime=(1/n)*(t(cbind(Xmat1,Xmat2,Xmat3))%*%cbind(Xmat1,Xmat2,Xmat3))
  
  
  n=dim(Xmat1)[1]
  p=ncol(Xmat3)
  p.oracle=ncol(Xmat3.oracle)
  
  bdy.knots.b.1 = c(0, max(y1 )) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2 -y1 ))
  
  b.1=BP.simData$b.1.bp
  b.2=BP.simData$b.2.bp
  b.3.y2my1=BP.simData$b.3.y2my1.bp
  #extracting starting values based on bernstein polynomials:
  startVals=BP.simData$startVals.bp
  startVals.oracle=startVals[nonzero.param.truth.index]
  
  
  fit.bp=FreqID.LT.bSpline.bp(Y=Y, lin.pred=lin.pred, data=data, startVals=startVals, frailty=TRUE, 
                              b.1, b.2, b.3.y2my1, 
                              bdy.knots.b.1, 
                              bdy.knots.b.2, 
                              bdy.knots.b.3.y2my1,
                              method)
  print("I successfully ran freqID.LT for the complete data, and now I am running it for finding oracle unpenalized estimate of oracle")
  if ('try-error' %in% class(fit.bp)){
    cat("I have a problem in running freqID.LT in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp=fit.bp
  }
  
  if (risks.status=="fullyshared"){
    fit.bp.oracle=FreqID.LT.bSpline.bp(Y=Y.oracle, lin.pred=lin.pred.oracle, 
                                       data=data.oracle.1,
                                       startVals=startVals.oracle , frailty,
                                       b.1, b.2, b.3.y2my1,
                                       bdy.knots.b.1,
                                       bdy.knots.b.2,
                                       bdy.knots.b.3.y2my1,
                                       method)
  }
  if (risks.status=="overlapping"){
    fit.bp.oracle=FreqID.LT.bSpline.bp.oracle.overlapping(Y=Y.oracle.1, lin.pred=lin.pred.oracle,
                                                          data.oracle.1,data.oracle.2,data.oracle.3,           
                                                          startVals=startVals.oracle , frailty,
                                                          b.1, b.2, b.3.y2my1,
                                                          bdy.knots.b.1,
                                                          bdy.knots.b.2,
                                                          bdy.knots.b.3.y2my1,
                                                          method)    
  }
  
  
  print("I successfully ran freqID.LT.oracle for the complete data, and results are coming out shortly")
  
  #If cannot fit a model to find the unpenalized estimate, just jump over this iteration, and go to the next!
  
  
  if ('try-error' %in% class(fit.bp.oracle)){
    cat("I have a problem in running freqID.LT.oracle in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp.oracle=fit.bp.oracle
  }
  
  unpen.est=(fit.bp$estimate)[(nuisance.p+1):(length(BP.simData$startVals.bp))]
  unpen.est.oracle=(fit.bp.oracle$estimate)[(nuisance.p+1):(length(nonzero.param.truth.index))]
  
  nuisance.estimates=(fit.bp$estimate)[1:nuisance.p]
  nuisance.estimates.oracle=(fit.bp.oracle$estimate)[1:nuisance.p]
  
  E_ZZprimei.oracle=list(0,n)
  
  for (i in 1:n){
    E_ZZprimei.oracle[[i]]=matrix(cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1)%*%t(matrix (cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1))
  }
  E_ZZprime.oracle=(1/n)*Reduce(`+`, E_ZZprimei.oracle)
  
  return(list(fit.bp=fit.bp,fit.bp.oracle.estimate=fit.bp.oracle$estimate,data=data,BP.simData=BP.simData,
              unpen.est.oracle=unpen.est.oracle, unpen.est=unpen.est,E_ZZprime=E_ZZprime,
              E_ZZprime.oracle=E_ZZprime.oracle,Y=Y,Y.oracle=Y.oracle,
              Y.oracle.1=Y.oracle.1,Y.oracle.2=Y.oracle.2,Y.oracle.3=Y.oracle.3,
              Xmat1.oracle=Xmat1.oracle,Xmat2.oracle=Xmat2.oracle,Xmat3.oracle=Xmat3.oracle,
              lin.pred=lin.pred,lin.pred.oracle=lin.pred.oracle,data=data,
              data.oracle.1=data.oracle.1,data.oracle.2=data.oracle.2,data.oracle.3=data.oracle.3,
              y1=y1,y2=y2,delta1=delta1,delta2=delta2,l=l,nuisance.estimates=nuisance.estimates,
              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,beta.truth.oracle=beta.truth.oracle,beta.truth=beta.truth
              ,startVals=startVals,startVals.oracle=startVals.oracle,b.1=b.1,b.2=b.2,b.3.y2my1=b.3.y2my1 ))
}

data.generator.and.unpenalized.estimator.12cov=function(turn,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution){
  set.seed(turn)
  # # beta1.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta2.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta3.true=c(-1,-1,-1,-1,0,0,0,0,0,0)
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  
  
  if(grp.or.indiv=="indiv"){
    BP.simData=genWB.simData.12cov.BP.BS(mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution)
  }
  
  
  p=(dim(BP.simData$Y)[2])-5
  #dimension when oracle:
  beta1.true.nonzero.length=sum(beta1.true!=0)
  beta2.true.nonzero.length=sum(beta2.true!=0)
  beta3.true.nonzero.length=sum(beta3.true!=0)
  p.oracle=sum(beta1.true.nonzero.length,beta2.true.nonzero.length,beta3.true.nonzero.length)
  
  #true beta values:
  beta.truth=(BP.simData$WBpara)[8:((3*p)+7)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$WBpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$WBpara)[8:((3*p)+7)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$WBpara)[nonzero.param.truth.index]
  ## Data and model specifications
  #true beta values:
  nuisance.p= sum(BP.simData$m)+3+1
  beta.truth=(BP.simData$BSpara)[(nuisance.p+1):((3*p)+nuisance.p)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$BSpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$BPpara)[(nuisance.p+1):((3*p)+nuisance.p)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$BPpara)[nonzero.param.truth.index]
  ## Data and model specifications
  Y <- BP.simData$Y                    # y1, delta1, y2, delta2, l = left-truncation time
  
  Y.oracle=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
  if (risks.status=="fullyshared"){
    Y.oracle.1=Y.oracle
    Y.oracle.2=Y.oracle
    Y.oracle.3=Y.oracle
  }
  if(risks.status=="overlapping"){
    Y.oracle.1=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
    Y.oracle.2=(BP.simData$Y)[,-(5+(which(beta2.true==0)))]
    Y.oracle.3=(BP.simData$Y)[,-(5+(which(beta3.true==0)))]
  }
  
  lin.pred <- BP.simData$lin.pred      # linear predictor formula
  lin.pred.oracle =BP.simData$lin.pred.oracle
  
  data <- BP.simData$data              # simulated data
  
  if (risks.status=="fullyshared"){
    data.oracle.1=data.oracle.2=data.oracle.3=BP.simData$data[,-(5+(which(beta1.true==0)))]
  }
  if (risks.status=="overlapping"){
    data.oracle.1=BP.simData$data[,-(5+(which(beta1.true==0)))]
    data.oracle.2=BP.simData$data[,-(5+(which(beta2.true==0)))]
    data.oracle.3=BP.simData$data[,-(5+(which(beta3.true==0)))]
    
  }
  
  trueBPparams <- BP.simData$BPpara    # true parameters used in data generation
  trueBPparams.oracle <- BP.simData$BPpara[nonzero.param.truth.index]
  frailty <- BP.simData$frailty        # allowing for a shared frailty
  para=BP.simData$BPpara
  para.oracle=trueBPparams.oracle 
  y1=BP.simData$Y[,1]
  y2=BP.simData$Y[,3]
  delta1=BP.simData$Y[,2]
  delta2=BP.simData$Y[,4]
  l=BP.simData$Y[,5]
  Ylength=dim(Y)[2]
  Ylength.oracle=dim(Y.oracle)[2]
  Y[,(6:Ylength)]->CovMat
  Y.oracle[,(6:Ylength.oracle)]->CovMat.oracle
  
  
  Xmat1=Xmat2=Xmat3=as.matrix(CovMat)
  if (risks.status=="fullyshared"){
    Xmat1.oracle=Xmat2.oracle=Xmat3.oracle=as.matrix(CovMat.oracle)
  }
  p=ncol(Xmat1)
  if(risks.status=="overlapping"){
    Xmat1.oracle=as.matrix(CovMat[,nonzero.beta.truth.index[1:4]])
    Xmat2.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[5:8])-p)])
    Xmat3.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[9:12])-(2*p))])
    
  }
  
  
  E_ZZprime=(1/n)*(t(cbind(Xmat1,Xmat2,Xmat3))%*%cbind(Xmat1,Xmat2,Xmat3))
  
  
  n=dim(Xmat1)[1]
  p=ncol(Xmat3)
  p.oracle=ncol(Xmat3.oracle)
  
  bdy.knots.b.1 = c(0, max(y1 )) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2 -y1 ))
  
  b.1=BP.simData$b.1.bp
  b.2=BP.simData$b.2.bp
  b.3.y2my1=BP.simData$b.3.y2my1.bp
  #extracting starting values based on bernstein polynomials:
  startVals=BP.simData$startVals.bp
  startVals.oracle=startVals[nonzero.param.truth.index]
  
  
  fit.bp=FreqID.LT.bSpline.bp(Y=Y, lin.pred=lin.pred, data=data, startVals=startVals, frailty=TRUE, 
                              b.1, b.2, b.3.y2my1, 
                              bdy.knots.b.1, 
                              bdy.knots.b.2, 
                              bdy.knots.b.3.y2my1,
                              method)
  print("I successfully ran freqID.LT for the complete data, and now I am running it for finding oracle unpenalized estimate of oracle")
  if ('try-error' %in% class(fit.bp)){
    cat("I have a problem in running freqID.LT in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp=fit.bp
  }
  
  if (risks.status=="fullyshared"){
    fit.bp.oracle=FreqID.LT.bSpline.bp(Y=Y.oracle, lin.pred=lin.pred.oracle, 
                                       data=data.oracle.1,
                                       startVals=startVals.oracle , frailty,
                                       b.1, b.2, b.3.y2my1,
                                       bdy.knots.b.1,
                                       bdy.knots.b.2,
                                       bdy.knots.b.3.y2my1,
                                       method)
  }
  if (risks.status=="overlapping"){
    fit.bp.oracle=FreqID.LT.bSpline.bp.oracle.overlapping(Y=Y.oracle.1, lin.pred=lin.pred.oracle,
                                                          data.oracle.1,data.oracle.2,data.oracle.3,           
                                                          startVals=startVals.oracle , frailty,
                                                          b.1, b.2, b.3.y2my1,
                                                          bdy.knots.b.1,
                                                          bdy.knots.b.2,
                                                          bdy.knots.b.3.y2my1,
                                                          method)    
  }
  
  
  print("I successfully ran freqID.LT.oracle for the complete data, and results are coming out shortly")
  
  #If cannot fit a model to find the unpenalized estimate, just jump over this iteration, and go to the next!
  
  
  if ('try-error' %in% class(fit.bp.oracle)){
    cat("I have a problem in running freqID.LT.oracle in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp.oracle=fit.bp.oracle
  }
  
  unpen.est=(fit.bp$estimate)[(nuisance.p+1):(length(BP.simData$startVals.bp))]
  unpen.est.oracle=(fit.bp.oracle$estimate)[(nuisance.p+1):(length(nonzero.param.truth.index))]
  
  nuisance.estimates=(fit.bp$estimate)[1:nuisance.p]
  nuisance.estimates.oracle=(fit.bp.oracle$estimate)[1:nuisance.p]
  
  E_ZZprimei.oracle=list(0,n)
  
  for (i in 1:n){
    E_ZZprimei.oracle[[i]]=matrix(cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1)%*%t(matrix (cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1))
  }
  E_ZZprime.oracle=(1/n)*Reduce(`+`, E_ZZprimei.oracle)
  
  return(list(fit.bp=fit.bp,fit.bp.oracle.estimate=fit.bp.oracle$estimate,data=data,BP.simData=BP.simData,
              unpen.est.oracle=unpen.est.oracle, unpen.est=unpen.est,E_ZZprime=E_ZZprime,
              E_ZZprime.oracle=E_ZZprime.oracle,Y=Y,Y.oracle=Y.oracle,
              Y.oracle.1=Y.oracle.1,Y.oracle.2=Y.oracle.2,Y.oracle.3=Y.oracle.3,
              Xmat1.oracle=Xmat1.oracle,Xmat2.oracle=Xmat2.oracle,Xmat3.oracle=Xmat3.oracle,
              lin.pred=lin.pred,lin.pred.oracle=lin.pred.oracle,data=data,
              data.oracle.1=data.oracle.1,data.oracle.2=data.oracle.2,data.oracle.3=data.oracle.3,
              y1=y1,y2=y2,delta1=delta1,delta2=delta2,l=l,nuisance.estimates=nuisance.estimates,
              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,beta.truth.oracle=beta.truth.oracle,beta.truth=beta.truth
              ,startVals=startVals,startVals.oracle=startVals.oracle,b.1=b.1,b.2=b.2,b.3.y2my1=b.3.y2my1))
}
data.generator.and.unpenalized.estimator.15cov=function(turn,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution){
  set.seed(turn)
  # # beta1.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta2.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta3.true=c(-1,-1,-1,-1,0,0,0,0,0,0)
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  
  
  if(grp.or.indiv=="indiv"){
    BP.simData=genWB.simData.15cov.BP.BS(mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution)
  }
  
  
  p=(dim(BP.simData$Y)[2])-5
  #dimension when oracle:
  beta1.true.nonzero.length=sum(beta1.true!=0)
  beta2.true.nonzero.length=sum(beta2.true!=0)
  beta3.true.nonzero.length=sum(beta3.true!=0)
  p.oracle=sum(beta1.true.nonzero.length,beta2.true.nonzero.length,beta3.true.nonzero.length)
  
  #true beta values:
  beta.truth=(BP.simData$WBpara)[8:((3*p)+7)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$WBpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$WBpara)[8:((3*p)+7)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$WBpara)[nonzero.param.truth.index]
  ## Data and model specifications
  #true beta values:
  nuisance.p= sum(BP.simData$m)+3+1
  beta.truth=(BP.simData$BSpara)[(nuisance.p+1):((3*p)+nuisance.p)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$BSpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$BPpara)[(nuisance.p+1):((3*p)+nuisance.p)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$BPpara)[nonzero.param.truth.index]
  ## Data and model specifications
  Y <- BP.simData$Y                    # y1, delta1, y2, delta2, l = left-truncation time
  
  Y.oracle=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
  if (risks.status=="fullyshared"){
    Y.oracle.1=Y.oracle
    Y.oracle.2=Y.oracle
    Y.oracle.3=Y.oracle
  }
  if(risks.status=="overlapping"){
    Y.oracle.1=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
    Y.oracle.2=(BP.simData$Y)[,-(5+(which(beta2.true==0)))]
    Y.oracle.3=(BP.simData$Y)[,-(5+(which(beta3.true==0)))]
  }
  
  lin.pred <- BP.simData$lin.pred      # linear predictor formula
  lin.pred.oracle =BP.simData$lin.pred.oracle
  
  data <- BP.simData$data              # simulated data
  
  if (risks.status=="fullyshared"){
    data.oracle.1=data.oracle.2=data.oracle.3=BP.simData$data[,-(5+(which(beta1.true==0)))]
  }
  if (risks.status=="overlapping"){
    data.oracle.1=BP.simData$data[,-(5+(which(beta1.true==0)))]
    data.oracle.2=BP.simData$data[,-(5+(which(beta2.true==0)))]
    data.oracle.3=BP.simData$data[,-(5+(which(beta3.true==0)))]
    
  }
  
  trueBPparams <- BP.simData$BPpara    # true parameters used in data generation
  trueBPparams.oracle <- BP.simData$BPpara[nonzero.param.truth.index]
  frailty <- BP.simData$frailty        # allowing for a shared frailty
  para=BP.simData$BPpara
  para.oracle=trueBPparams.oracle 
  y1=BP.simData$Y[,1]
  y2=BP.simData$Y[,3]
  delta1=BP.simData$Y[,2]
  delta2=BP.simData$Y[,4]
  l=BP.simData$Y[,5]
  Ylength=dim(Y)[2]
  Ylength.oracle=dim(Y.oracle)[2]
  Y[,(6:Ylength)]->CovMat
  Y.oracle[,(6:Ylength.oracle)]->CovMat.oracle
  
  
  Xmat1=Xmat2=Xmat3=as.matrix(CovMat)
  if (risks.status=="fullyshared"){
    Xmat1.oracle=Xmat2.oracle=Xmat3.oracle=as.matrix(CovMat.oracle)
  }
  p=ncol(Xmat1)
  if(risks.status=="overlapping"){
    Xmat1.oracle=as.matrix(CovMat[,nonzero.beta.truth.index[1:4]])
    Xmat2.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[5:8])-p)])
    Xmat3.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[9:12])-(2*p))])
    
  }
  
  
  E_ZZprime=(1/n)*(t(cbind(Xmat1,Xmat2,Xmat3))%*%cbind(Xmat1,Xmat2,Xmat3))
  
  
  n=dim(Xmat1)[1]
  p=ncol(Xmat3)
  p.oracle=ncol(Xmat3.oracle)
  
  bdy.knots.b.1 = c(0, max(y1 )) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2 -y1 ))
  
  b.1=BP.simData$b.1.bp
  b.2=BP.simData$b.2.bp
  b.3.y2my1=BP.simData$b.3.y2my1.bp
  #extracting starting values based on bernstein polynomials:
  startVals=BP.simData$startVals.bp
  startVals.oracle=startVals[nonzero.param.truth.index]
  
  
  fit.bp=FreqID.LT.bSpline.bp(Y=Y, lin.pred=lin.pred, data=data, startVals=startVals, frailty=TRUE, 
                              b.1, b.2, b.3.y2my1, 
                              bdy.knots.b.1, 
                              bdy.knots.b.2, 
                              bdy.knots.b.3.y2my1,
                              method)
  print("I successfully ran freqID.LT for the complete data, and now I am running it for finding oracle unpenalized estimate of oracle")
  if ('try-error' %in% class(fit.bp)){
    cat("I have a problem in running freqID.LT in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp=fit.bp
  }
  
  if (risks.status=="fullyshared"){
    fit.bp.oracle=FreqID.LT.bSpline.bp(Y=Y.oracle, lin.pred=lin.pred.oracle, 
                                       data=data.oracle.1,
                                       startVals=startVals.oracle , frailty,
                                       b.1, b.2, b.3.y2my1,
                                       bdy.knots.b.1,
                                       bdy.knots.b.2,
                                       bdy.knots.b.3.y2my1,
                                       method)
  }
  if (risks.status=="overlapping"){
    fit.bp.oracle=FreqID.LT.bSpline.bp.oracle.overlapping(Y=Y.oracle.1, lin.pred=lin.pred.oracle,
                                                          data.oracle.1,data.oracle.2,data.oracle.3,           
                                                          startVals=startVals.oracle , frailty,
                                                          b.1, b.2, b.3.y2my1,
                                                          bdy.knots.b.1,
                                                          bdy.knots.b.2,
                                                          bdy.knots.b.3.y2my1,
                                                          method)    
  }
  
  
  print("I successfully ran freqID.LT.oracle for the complete data, and results are coming out shortly")
  
  #If cannot fit a model to find the unpenalized estimate, just jump over this iteration, and go to the next!
  
  
  if ('try-error' %in% class(fit.bp.oracle)){
    cat("I have a problem in running freqID.LT.oracle in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp.oracle=fit.bp.oracle
  }
  
  unpen.est=(fit.bp$estimate)[(nuisance.p+1):(length(BP.simData$startVals.bp))]
  unpen.est.oracle=(fit.bp.oracle$estimate)[(nuisance.p+1):(length(nonzero.param.truth.index))]
  
  nuisance.estimates=(fit.bp$estimate)[1:nuisance.p]
  nuisance.estimates.oracle=(fit.bp.oracle$estimate)[1:nuisance.p]
  
  E_ZZprimei.oracle=list(0,n)
  
  for (i in 1:n){
    E_ZZprimei.oracle[[i]]=matrix(cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1)%*%t(matrix (cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1))
  }
  E_ZZprime.oracle=(1/n)*Reduce(`+`, E_ZZprimei.oracle)
  
  return(list(fit.bp=fit.bp,fit.bp.oracle.estimate=fit.bp.oracle$estimate,data=data,BP.simData=BP.simData,
              unpen.est.oracle=unpen.est.oracle, unpen.est=unpen.est,E_ZZprime=E_ZZprime,
              E_ZZprime.oracle=E_ZZprime.oracle,Y=Y,Y.oracle=Y.oracle,
              Y.oracle.1=Y.oracle.1,Y.oracle.2=Y.oracle.2,Y.oracle.3=Y.oracle.3,
              Xmat1.oracle=Xmat1.oracle,Xmat2.oracle=Xmat2.oracle,Xmat3.oracle=Xmat3.oracle,
              lin.pred=lin.pred,lin.pred.oracle=lin.pred.oracle,data=data,
              data.oracle.1=data.oracle.1,data.oracle.2=data.oracle.2,data.oracle.3=data.oracle.3,
              y1=y1,y2=y2,delta1=delta1,delta2=delta2,l=l,nuisance.estimates=nuisance.estimates,
              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,beta.truth.oracle=beta.truth.oracle,beta.truth=beta.truth
              ,startVals=startVals,startVals.oracle=startVals.oracle,b.1=b.1,b.2=b.2,b.3.y2my1=b.3.y2my1))
}
data.generator.and.unpenalized.estimator.16cov=function(turn,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution){
  set.seed(turn)
  # # beta1.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta2.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta3.true=c(-1,-1,-1,-1,0,0,0,0,0,0)
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  
  
  if(grp.or.indiv=="indiv"){
    BP.simData=genWB.simData.16cov.BP.BS(mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution)
  }
  
  
  p=(dim(BP.simData$Y)[2])-5
  #dimension when oracle:
  beta1.true.nonzero.length=sum(beta1.true!=0)
  beta2.true.nonzero.length=sum(beta2.true!=0)
  beta3.true.nonzero.length=sum(beta3.true!=0)
  p.oracle=sum(beta1.true.nonzero.length,beta2.true.nonzero.length,beta3.true.nonzero.length)
  
  #true beta values:
  beta.truth=(BP.simData$WBpara)[8:((3*p)+7)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$WBpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$WBpara)[8:((3*p)+7)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$WBpara)[nonzero.param.truth.index]
  ## Data and model specifications
  #true beta values:
  nuisance.p= sum(BP.simData$m)+3+1
  beta.truth=(BP.simData$BSpara)[(nuisance.p+1):((3*p)+nuisance.p)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$BSpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$BPpara)[(nuisance.p+1):((3*p)+nuisance.p)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$BPpara)[nonzero.param.truth.index]
  ## Data and model specifications
  Y <- BP.simData$Y                    # y1, delta1, y2, delta2, l = left-truncation time
  
  Y.oracle=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
  if (risks.status=="fullyshared"){
    Y.oracle.1=Y.oracle
    Y.oracle.2=Y.oracle
    Y.oracle.3=Y.oracle
  }
  if(risks.status=="overlapping"){
    Y.oracle.1=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
    Y.oracle.2=(BP.simData$Y)[,-(5+(which(beta2.true==0)))]
    Y.oracle.3=(BP.simData$Y)[,-(5+(which(beta3.true==0)))]
  }
  
  lin.pred <- BP.simData$lin.pred      # linear predictor formula
  lin.pred.oracle =BP.simData$lin.pred.oracle
  
  data <- BP.simData$data              # simulated data
  
  if (risks.status=="fullyshared"){
    data.oracle.1=data.oracle.2=data.oracle.3=BP.simData$data[,-(5+(which(beta1.true==0)))]
  }
  if (risks.status=="overlapping"){
    data.oracle.1=BP.simData$data[,-(5+(which(beta1.true==0)))]
    data.oracle.2=BP.simData$data[,-(5+(which(beta2.true==0)))]
    data.oracle.3=BP.simData$data[,-(5+(which(beta3.true==0)))]
    
  }
  
  trueBPparams <- BP.simData$BPpara    # true parameters used in data generation
  trueBPparams.oracle <- BP.simData$BPpara[nonzero.param.truth.index]
  frailty <- BP.simData$frailty        # allowing for a shared frailty
  para=BP.simData$BPpara
  para.oracle=trueBPparams.oracle 
  y1=BP.simData$Y[,1]
  y2=BP.simData$Y[,3]
  delta1=BP.simData$Y[,2]
  delta2=BP.simData$Y[,4]
  l=BP.simData$Y[,5]
  Ylength=dim(Y)[2]
  Ylength.oracle=dim(Y.oracle)[2]
  Y[,(6:Ylength)]->CovMat
  Y.oracle[,(6:Ylength.oracle)]->CovMat.oracle
  
  
  Xmat1=Xmat2=Xmat3=as.matrix(CovMat)
  if (risks.status=="fullyshared"){
    Xmat1.oracle=Xmat2.oracle=Xmat3.oracle=as.matrix(CovMat.oracle)
  }
  p=ncol(Xmat1)
  if(risks.status=="overlapping"){
    Xmat1.oracle=as.matrix(CovMat[,nonzero.beta.truth.index[1:4]])
    Xmat2.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[5:8])-p)])
    Xmat3.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[9:12])-(2*p))])
    
  }
  
  
  E_ZZprime=(1/n)*(t(cbind(Xmat1,Xmat2,Xmat3))%*%cbind(Xmat1,Xmat2,Xmat3))
  
  
  n=dim(Xmat1)[1]
  p=ncol(Xmat3)
  p.oracle=ncol(Xmat3.oracle)
  
  bdy.knots.b.1 = c(0, max(y1 )) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2 -y1 ))
  
  b.1=BP.simData$b.1.bp
  b.2=BP.simData$b.2.bp
  b.3.y2my1=BP.simData$b.3.y2my1.bp
  #extracting starting values based on bernstein polynomials:
  startVals=BP.simData$startVals.bp
  startVals.oracle=startVals[nonzero.param.truth.index]
  
  
  fit.bp=FreqID.LT.bSpline.bp(Y=Y, lin.pred=lin.pred, data=data, startVals=startVals, frailty=TRUE, 
                              b.1, b.2, b.3.y2my1, 
                              bdy.knots.b.1, 
                              bdy.knots.b.2, 
                              bdy.knots.b.3.y2my1,
                              method)
  print("I successfully ran freqID.LT for the complete data, and now I am running it for finding oracle unpenalized estimate of oracle")
  if ('try-error' %in% class(fit.bp)){
    cat("I have a problem in running freqID.LT in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp=fit.bp
  }
  
  if (risks.status=="fullyshared"){
    fit.bp.oracle=FreqID.LT.bSpline.bp(Y=Y.oracle, lin.pred=lin.pred.oracle, 
                                       data=data.oracle.1,
                                       startVals=startVals.oracle , frailty,
                                       b.1, b.2, b.3.y2my1,
                                       bdy.knots.b.1,
                                       bdy.knots.b.2,
                                       bdy.knots.b.3.y2my1,
                                       method)
  }
  if (risks.status=="overlapping"){
    fit.bp.oracle=FreqID.LT.bSpline.bp.oracle.overlapping(Y=Y.oracle.1, lin.pred=lin.pred.oracle,
                                                          data.oracle.1,data.oracle.2,data.oracle.3,           
                                                          startVals=startVals.oracle , frailty,
                                                          b.1, b.2, b.3.y2my1,
                                                          bdy.knots.b.1,
                                                          bdy.knots.b.2,
                                                          bdy.knots.b.3.y2my1,
                                                          method)    
  }
  
  
  print("I successfully ran freqID.LT.oracle for the complete data, and results are coming out shortly")
  
  #If cannot fit a model to find the unpenalized estimate, just jump over this iteration, and go to the next!
  
  
  if ('try-error' %in% class(fit.bp.oracle)){
    cat("I have a problem in running freqID.LT.oracle in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp.oracle=fit.bp.oracle
  }
  
  unpen.est=(fit.bp$estimate)[(nuisance.p+1):(length(BP.simData$startVals.bp))]
  unpen.est.oracle=(fit.bp.oracle$estimate)[(nuisance.p+1):(length(nonzero.param.truth.index))]
  
  nuisance.estimates=(fit.bp$estimate)[1:nuisance.p]
  nuisance.estimates.oracle=(fit.bp.oracle$estimate)[1:nuisance.p]
  
  E_ZZprimei.oracle=list(0,n)
  
  for (i in 1:n){
    E_ZZprimei.oracle[[i]]=matrix(cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1)%*%t(matrix (cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1))
  }
  E_ZZprime.oracle=(1/n)*Reduce(`+`, E_ZZprimei.oracle)
  
  return(list(fit.bp=fit.bp,fit.bp.oracle.estimate=fit.bp.oracle$estimate,data=data,BP.simData=BP.simData,
              unpen.est.oracle=unpen.est.oracle, unpen.est=unpen.est,E_ZZprime=E_ZZprime,
              E_ZZprime.oracle=E_ZZprime.oracle,Y=Y,Y.oracle=Y.oracle,
              Y.oracle.1=Y.oracle.1,Y.oracle.2=Y.oracle.2,Y.oracle.3=Y.oracle.3,
              Xmat1.oracle=Xmat1.oracle,Xmat2.oracle=Xmat2.oracle,Xmat3.oracle=Xmat3.oracle,
              lin.pred=lin.pred,lin.pred.oracle=lin.pred.oracle,data=data,
              data.oracle.1=data.oracle.1,data.oracle.2=data.oracle.2,data.oracle.3=data.oracle.3,
              y1=y1,y2=y2,delta1=delta1,delta2=delta2,l=l,nuisance.estimates=nuisance.estimates,
              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,beta.truth.oracle=beta.truth.oracle,beta.truth=beta.truth
              ,startVals=startVals,startVals.oracle=startVals.oracle,b.1=b.1,b.2=b.2,b.3.y2my1=b.3.y2my1))
}

data.generator.and.unpenalized.estimator.20cov=function(turn,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution){
  set.seed(turn)
  # # beta1.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta2.true=c(1,1,1,1,0,0,0,0,0,0)
  # # beta3.true=c(-1,-1,-1,-1,0,0,0,0,0,0)
  # beta1.true = c(-0.8,1,1,0.9,0,0,0,0,0,0)
  # beta2.true = c(1,1,1,0.9,0,0,0,0,0,0)
  # beta3.true = c(-1,1,0.9,1,0,0,0,0,0,0)
  
  
  if(grp.or.indiv=="indiv"){
    BP.simData=genWB.simData.20cov.BP.BS(mu1,mu2,sigma1,sigma2,p1,p2,weibull.param.log,c1,c2,lt1,lt2,n ,m,rho,beta1.true,beta2.true,beta3.true, frailty , risks.status,yvalslength,theta.given,frailty_distribution)
  }
  
  
  p=(dim(BP.simData$Y)[2])-5
  #dimension when oracle:
  beta1.true.nonzero.length=sum(beta1.true!=0)
  beta2.true.nonzero.length=sum(beta2.true!=0)
  beta3.true.nonzero.length=sum(beta3.true!=0)
  p.oracle=sum(beta1.true.nonzero.length,beta2.true.nonzero.length,beta3.true.nonzero.length)
  
  #true beta values:
  beta.truth=(BP.simData$WBpara)[8:((3*p)+7)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$WBpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$WBpara)[8:((3*p)+7)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$WBpara)[nonzero.param.truth.index]
  ## Data and model specifications
  #true beta values:
  nuisance.p= sum(BP.simData$m)+3+1
  beta.truth=(BP.simData$BSpara)[(nuisance.p+1):((3*p)+nuisance.p)]
  #oracle true beta values:
  nonzero.param.truth.index=which((BP.simData$BSpara)!=0)
  nonzero.beta.truth.index=which(beta.truth!=0)
  beta.truth.oracle=((BP.simData$BPpara)[(nuisance.p+1):((3*p)+nuisance.p)])[nonzero.beta.truth.index]
  param.truth.oracle=(BP.simData$BPpara)[nonzero.param.truth.index]
  ## Data and model specifications
  Y <- BP.simData$Y                    # y1, delta1, y2, delta2, l = left-truncation time
  
  Y.oracle=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
  if (risks.status=="fullyshared"){
    Y.oracle.1=Y.oracle
    Y.oracle.2=Y.oracle
    Y.oracle.3=Y.oracle
  }
  if(risks.status=="overlapping"){
    Y.oracle.1=(BP.simData$Y)[,-(5+(which(beta1.true==0)))]
    Y.oracle.2=(BP.simData$Y)[,-(5+(which(beta2.true==0)))]
    Y.oracle.3=(BP.simData$Y)[,-(5+(which(beta3.true==0)))]
  }
  
  lin.pred <- BP.simData$lin.pred      # linear predictor formula
  lin.pred.oracle =BP.simData$lin.pred.oracle
  
  data <- BP.simData$data              # simulated data
  
  if (risks.status=="fullyshared"){
    data.oracle.1=data.oracle.2=data.oracle.3=BP.simData$data[,-(5+(which(beta1.true==0)))]
  }
  if (risks.status=="overlapping"){
    data.oracle.1=BP.simData$data[,-(5+(which(beta1.true==0)))]
    data.oracle.2=BP.simData$data[,-(5+(which(beta2.true==0)))]
    data.oracle.3=BP.simData$data[,-(5+(which(beta3.true==0)))]
    
  }
  
  trueBPparams <- BP.simData$BPpara    # true parameters used in data generation
  trueBPparams.oracle <- BP.simData$BPpara[nonzero.param.truth.index]
  frailty <- BP.simData$frailty        # allowing for a shared frailty
  para=BP.simData$BPpara
  para.oracle=trueBPparams.oracle 
  y1=BP.simData$Y[,1]
  y2=BP.simData$Y[,3]
  delta1=BP.simData$Y[,2]
  delta2=BP.simData$Y[,4]
  l=BP.simData$Y[,5]
  Ylength=dim(Y)[2]
  Ylength.oracle=dim(Y.oracle)[2]
  Y[,(6:Ylength)]->CovMat
  Y.oracle[,(6:Ylength.oracle)]->CovMat.oracle
  
  
  Xmat1=Xmat2=Xmat3=as.matrix(CovMat)
  if (risks.status=="fullyshared"){
    Xmat1.oracle=Xmat2.oracle=Xmat3.oracle=as.matrix(CovMat.oracle)
  }
  p=ncol(Xmat1)
  if(risks.status=="overlapping"){
    Xmat1.oracle=as.matrix(CovMat[,nonzero.beta.truth.index[1:4]])
    Xmat2.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[5:8])-p)])
    Xmat3.oracle=as.matrix(CovMat[,((nonzero.beta.truth.index[9:12])-(2*p))])
    
  }
  
  
  E_ZZprime=(1/n)*(t(cbind(Xmat1,Xmat2,Xmat3))%*%cbind(Xmat1,Xmat2,Xmat3))
  
  
  n=dim(Xmat1)[1]
  p=ncol(Xmat3)
  p.oracle=ncol(Xmat3.oracle)
  
  bdy.knots.b.1 = c(0, max(y1 )) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2 -y1 ))
  
  b.1=BP.simData$b.1.bp
  b.2=BP.simData$b.2.bp
  b.3.y2my1=BP.simData$b.3.y2my1.bp
  #extracting starting values based on bernstein polynomials:
  startVals=BP.simData$startVals.bp
  startVals.oracle=startVals[nonzero.param.truth.index]
  
  
  fit.bp=FreqID.LT.bSpline.bp(Y=Y, lin.pred=lin.pred, data=data, startVals=startVals, frailty=TRUE, 
                              b.1, b.2, b.3.y2my1, 
                              bdy.knots.b.1, 
                              bdy.knots.b.2, 
                              bdy.knots.b.3.y2my1,
                              method)
  print("I successfully ran freqID.LT for the complete data, and now I am running it for finding oracle unpenalized estimate of oracle")
  if ('try-error' %in% class(fit.bp)){
    cat("I have a problem in running freqID.LT in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp=fit.bp
  }
  
  if (risks.status=="fullyshared"){
    fit.bp.oracle=FreqID.LT.bSpline.bp(Y=Y.oracle, lin.pred=lin.pred.oracle, 
                                       data=data.oracle.1,
                                       startVals=startVals.oracle , frailty,
                                       b.1, b.2, b.3.y2my1,
                                       bdy.knots.b.1,
                                       bdy.knots.b.2,
                                       bdy.knots.b.3.y2my1,
                                       method)
  }
  if (risks.status=="overlapping"){
    fit.bp.oracle=FreqID.LT.bSpline.bp.oracle.overlapping(Y=Y.oracle.1, lin.pred=lin.pred.oracle,
                                                          data.oracle.1,data.oracle.2,data.oracle.3,           
                                                          startVals=startVals.oracle , frailty,
                                                          b.1, b.2, b.3.y2my1,
                                                          bdy.knots.b.1,
                                                          bdy.knots.b.2,
                                                          bdy.knots.b.3.y2my1,
                                                          method)    
  }
  
  
  print("I successfully ran freqID.LT.oracle for the complete data, and results are coming out shortly")
  
  #If cannot fit a model to find the unpenalized estimate, just jump over this iteration, and go to the next!
  
  
  if ('try-error' %in% class(fit.bp.oracle)){
    cat("I have a problem in running freqID.LT.oracle in data.generator function at seed:",turn,"PLEASE CHECK THIS ITERATIon \n")
    next
  }else{
    fit.bp.oracle=fit.bp.oracle
  }
  
  unpen.est=(fit.bp$estimate)[(nuisance.p+1):(length(BP.simData$startVals.bp))]
  unpen.est.oracle=(fit.bp.oracle$estimate)[(nuisance.p+1):(length(nonzero.param.truth.index))]
  
  nuisance.estimates=(fit.bp$estimate)[1:nuisance.p]
  nuisance.estimates.oracle=(fit.bp.oracle$estimate)[1:nuisance.p]
  
  E_ZZprimei.oracle=list(0,n)
  
  for (i in 1:n){
    E_ZZprimei.oracle[[i]]=matrix(cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1)%*%t(matrix (cbind(Xmat1.oracle,Xmat2.oracle,Xmat3.oracle)[i,],(3*p.oracle),1))
  }
  E_ZZprime.oracle=(1/n)*Reduce(`+`, E_ZZprimei.oracle)
  
  return(list(fit.bp=fit.bp,fit.bp.oracle.estimate=fit.bp.oracle$estimate,data=data,BP.simData=BP.simData,
              unpen.est.oracle=unpen.est.oracle, unpen.est=unpen.est,E_ZZprime=E_ZZprime,
              E_ZZprime.oracle=E_ZZprime.oracle,Y=Y,Y.oracle=Y.oracle,
              Y.oracle.1=Y.oracle.1,Y.oracle.2=Y.oracle.2,Y.oracle.3=Y.oracle.3,
              Xmat1.oracle=Xmat1.oracle,Xmat2.oracle=Xmat2.oracle,Xmat3.oracle=Xmat3.oracle,
              lin.pred=lin.pred,lin.pred.oracle=lin.pred.oracle,data=data,
              data.oracle.1=data.oracle.1,data.oracle.2=data.oracle.2,data.oracle.3=data.oracle.3,
              y1=y1,y2=y2,delta1=delta1,delta2=delta2,l=l,nuisance.estimates=nuisance.estimates,
              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,beta.truth.oracle=beta.truth.oracle,beta.truth=beta.truth
              ,startVals=startVals,startVals.oracle=startVals.oracle,b.1=b.1,b.2=b.2,b.3.y2my1=b.3.y2my1))
}



ss2 <- function(j,tmpb,Q,B)
{
  a <- sum(tmpb*Q[,j])-tmpb[j]*Q[j,j]
  s <- 2*(a-B[j])
  return(s)
}

#initi could be the unpenalized estimate: unpen.est
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


bar.finder.iter.solvebar.bp=function(tol,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
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
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    bar <-  solveBAR(Y_irls,X_irls,lam,xi)
    
    # lasso <-  solveLasso(Y_irls,X_irls,lam)
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
  return(list(para.est=para.est,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}

bar.finder.iter.solvebar.bp.grp=function(G1,G2,G3,G4,tol,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,
                                         Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data,
                                         model = "semi-Markov", frailty, startVals,lam,xi,
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
  
  G1.extended=c(G1,10+G1,20+G1)
  G2.extended=c(G2,10+G2,20+G2)
  G3.extended=c(G3,10+G3,20+G3)
  G4.extended=c(G4,10+G4,20+G4)
  
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


AdaLasso.finder.iter.solveAdaLasso.bp=function(nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,
                                               delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, 
                                               model = "semi-Markov", frailty=frailty, startVals,lam,
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
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    AdaLasso <-  solveAdaLasso(3*p,X_irls,Y_irls,unpen.est,1/abs(unpen.est),lam)
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
  return(list(para.est=para.est,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}
AdaLasso.finder.iter.solveAdaLasso.bp.grp=function(G1,G2,G3,G4,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,
                                                   y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, 
                                                   model = "semi-Markov", frailty=frailty, startVals,lam,
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
    # cat("para.est",para.est,"in iter:",count,"\n")
    X_irls=chol(H)
    Y_irls=forwardsolve(t(X_irls),H%*%betaold-G)
    # lamvec=seq(0,10,0.05)#the values for penalizing bar with. 
    AdaLasso <-  solveAdaLasso(3*p,X_irls,Y_irls,unpen.est,1/abs(unpen.est),lam)
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
  G1.extended=c(G1,10+G1,20+G1)
  G2.extended=c(G2,10+G2,20+G2)
  G3.extended=c(G3,10+G3,20+G3)
  G4.extended=c(G4,10+G4,20+G4)
  
  G1.perc=ifelse(sum(abs(beta[G1.extended])!=0)==length(G1.extended),1,0)
  G2.perc=ifelse(sum(abs(beta[G2.extended])!=0)==length(G2.extended),1,0)
  G3.perc=ifelse(sum(abs(beta[G3.extended])==0)==length(G3.extended),1,0)
  G4.perc=ifelse(sum(abs(beta[G4.extended])==0)==length(G4.extended),1,0)
  FN=length(which(abs(beta[nonzero.truth.index])==0))
  return(list(G1.perc=G1.perc,G2.perc=G2.perc,G3.perc=G3.perc,G4.perc=G4.perc,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}

lasso.finder.iter.solvelasso.bp=function(nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,
                                         delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
                                         frailty=frailty, startVals,lam,
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
  return(list(para.est=para.est,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}

lasso.finder.iter.solvelasso.bp.grp=function(G1,G2,G3,G4,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,
                                             y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data,
                                             model = "semi-Markov", frailty=frailty, startVals,lam,
                                             b.1,  
                                             b.2,  
                                             b.3.y2my1){
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
  
  G1.extended=c(G1,10+G1,20+G1)
  G2.extended=c(G2,10+G2,20+G2)
  G3.extended=c(G3,10+G3,20+G3)
  G4.extended=c(G4,10+G4,20+G4)
  
  G1.perc=ifelse(sum(abs(beta[G1.extended])!=0)==length(G1.extended),1,0)
  G2.perc=ifelse(sum(abs(beta[G2.extended])!=0)==length(G2.extended),1,0)
  G3.perc=ifelse(sum(abs(beta[G3.extended])==0)==length(G3.extended),1,0)
  G4.perc=ifelse(sum(abs(beta[G4.extended])==0)==length(G4.extended),1,0)
  return(list(G1.perc=G1.perc,G2.perc=G2.perc,G3.perc=G3.perc,G4.perc=G4.perc,TP=TP,FP=FP,FN=FN,beta=beta,H=H,G=G,X=X_irls,y=Y_irls,MSE=MSE,count=count))
  
}



GCV.finder.lasso.bp=function(lambdavec,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,
                             delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov",
                             frailty, startVals,
                             b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    lasso.est=lasso.finder.iter.solvelasso.bp(nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,
                                              delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
                                              frailty, startVals,lam,
                                              b.1,b.2,b.3.y2my1)
    TP.for.that.lam[which(lam==lambdavec)]=lasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=lasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=lasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=lasso.est$FN
    theta.for.that.lam.lasso[,(which(lam==lambdavec))]=lasso.est$para.est[1:11]
    
    
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
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  theta.final=theta.for.that.lam.lasso[,optimal.gcv]
  
  return(list(theta.final=theta.final,GCV.lasso=GCV.lasso,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.lasso,betavarsel.lasso=betavarsel.lasso))
}

GCV.finder.lasso.bp.grp=function(G1,G2,G3,G4,lambdavec,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,
                                 y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov",
                                 frailty, startVals,
                                 b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    lasso.est=lasso.finder.iter.solvelasso.bp.grp(G1,G2,G3,G4,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,
                                                  y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data,
                                                  model = "semi-Markov", frailty, startVals,lam,
                                                  b.1,  
                                                  b.2,  
                                                  b.3.y2my1)
    TP.for.that.lam[which(lam==lambdavec)]=lasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=lasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=lasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=lasso.est$FN
    G1.for.that.lam[which(lam==lambdavec)]=lasso.est$G1.perc
    G2.for.that.lam[which(lam==lambdavec)]=lasso.est$G2.perc
    G3.for.that.lam[which(lam==lambdavec)]=lasso.est$G3.perc
    G4.for.that.lam[which(lam==lambdavec)]=lasso.est$G4.perc
    
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
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  G1.extended=c(G1,10+G1,20+G1)
  G2.extended=c(G2,10+G2,20+G2)
  G3.extended=c(G3,10+G3,20+G3)
  G4.extended=c(G4,10+G4,20+G4)
  
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

GCV.finder.AdaLasso.bp=function(lambdavec,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,
                                delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov", frailty,
                                startVals,
                                b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    AdaLasso.est=AdaLasso.finder.iter.solveAdaLasso.bp(nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,
                                                       delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, 
                                                       model = "semi-Markov", frailty=frailty, startVals,lam,
                                                       b.1,b.2,b.3.y2my1)
    
    TP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FN
    theta.for.that.lam.AdaLasso[,(which(lam==lambdavec))]=AdaLasso.est$para.est[1:11]
    
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
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  theta.final=theta.for.that.lam.AdaLasso[,optimal.gcv]
  
  
  return(list(theta.final=theta.final,GCV.AdaLasso=GCV.AdaLasso,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.AdaLasso,betavarsel.AdaLasso=betavarsel.AdaLasso))
}
GCV.finder.AdaLasso.bp.grp=function(G1,G2,G3,G4,lambdavec,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,
                                    y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov",
                                    frailty, startVals,
                                    b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    AdaLasso.est=AdaLasso.finder.iter.solveAdaLasso.bp.grp(G1,G2,G3,G4,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,
                                                           y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, 
                                                           model = "semi-Markov", frailty=frailty, startVals,lam,
                                                           b.1,b.2,b.3.y2my1)
    TP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$FN
    
    G1.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G1.perc
    G2.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G2.perc
    G3.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G3.perc
    G4.for.that.lam[which(lam==lambdavec)]=AdaLasso.est$G4.perc
    
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
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  G1.extended=c(G1,10+G1,20+G1)
  G2.extended=c(G2,10+G2,20+G2)
  G3.extended=c(G3,10+G3,20+G3)
  G4.extended=c(G4,10+G4,20+G4)
  
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

GCV.finder.BAR.bp=function(tol,lambdavec,xi,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,
                           delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov",
                           frailty, startVals,
                           b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    cat("this lam running:",lam,"\n")
    bar.est=bar.finder.iter.solvebar.bp(tol,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data, model = "semi-Markov",
                                        frailty=frailty,lam,xi,b.1,b.2,b.3.y2my1)
    TP.for.that.lam[which(lam==lambdavec)]=bar.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=bar.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=bar.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=bar.est$FN
    theta.for.that.lam.bar[,(which(lam==lambdavec))]=bar.est$para.est[1:11]
    
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
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  theta.final=theta.for.that.lam.bar[,optimal.gcv]
  
  return(list(theta.final=theta.final,GCV.bar=GCV.bar,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.bar))
}

GCV.finder.BAR.bp.grp=function(G1,G2,G3,G4,tol,lambdavec,xi,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,
                               Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3,lin.pred,data, model = "semi-Markov",
                               frailty, startVals,
                               b.1,b.2,b.3.y2my1){
  for (lam in lambdavec){
    bar.est=bar.finder.iter.solvebar.bp.grp(G1,G2,G3,G4,tol,nuisance.estimates,beta.truth,E_ZZprime,unpen.est,
                                            Y,y1,y2,delta1,delta2,l,Xmat1,Xmat2,Xmat3, lin.pred, data,
                                            model = "semi-Markov", frailty=frailty, startVals,lam,xi,
                                            b.1,b.2,b.3.y2my1)
    TP.for.that.lam[which(lam==lambdavec)]=bar.est$TP
    FP.for.that.lam[which(lam==lambdavec)]=bar.est$FP
    MSE.for.that.lam[which(lam==lambdavec)]=bar.est$MSE
    FN.for.that.lam[which(lam==lambdavec)]=bar.est$FN
    
    G1.for.that.lam[which(lam==lambdavec)]=bar.est$G1.perc
    G2.for.that.lam[which(lam==lambdavec)]=bar.est$G2.perc
    G3.for.that.lam[which(lam==lambdavec)]=bar.est$G3.perc
    G4.for.that.lam[which(lam==lambdavec)]=bar.est$G4.perc
    
    
    
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
  TP.final=TP.for.that.lam[optimal.gcv]
  FP.final=FP.for.that.lam[optimal.gcv]
  MSE.final=MSE.for.that.lam[optimal.gcv]
  FN.final=FN.for.that.lam[optimal.gcv]
  
  G1.extended=c(G1,10+G1,20+G1)
  G2.extended=c(G2,10+G2,20+G2)
  G3.extended=c(G3,10+G3,20+G3)
  G4.extended=c(G4,10+G4,20+G4)
  
  G1.perc.final=ifelse(sum(abs(beta.GCV[G1.extended])>tol)==length(G1.extended),1,0)
  G2.perc.final=ifelse(sum(abs(beta.GCV[G2.extended])>tol)==length(G2.extended),1,0)
  G3.perc.final=ifelse(sum(abs(beta.GCV[G3.extended])<=tol)==length(G3.extended),1,0)
  G4.perc.final=ifelse(sum(abs(beta.GCV[G4.extended])<=tol)==length(G4.extended),1,0)
  
  
  return(list(G1.perc.final=G1.perc.final,G2.perc.final=G2.perc.final,G3.perc.final=G3.perc.final,G4.perc.final=G4.perc.final,
              GCV.bar=GCV.bar,lam.final=lam.final,TP.final=TP.final,FP.final=FP.final,FN.final=FN.final,MSE.final=MSE.final,optimal.gcv=optimal.gcv,beta.GCV=beta.GCV,GCV=GCV.bar))
}

repeat.get.last.result.bp.update.theta=function(seeds,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,lambdavec.bar,lambdavec.lasso,lambdavec.ada,theta.given,frailty_distribution){
  for (f in seeds){
    cat("this is risks.stat:",risks.status,"\n")
    cat("Now running seed",f,"\n")
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method)
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status)
    if (length(beta1.true)==12){
      dataa=try(data.generator.and.unpenalized.estimator.12cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if (length(beta1.true)==15){
      dataa=try(data.generator.and.unpenalized.estimator.15cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if (length(beta1.true)==16){
      dataa=try(data.generator.and.unpenalized.estimator.16cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    
    if (length(beta1.true)==20){
      dataa=try(data.generator.and.unpenalized.estimator.20cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if (length(beta1.true)==10){
      dataa=try(data.generator.and.unpenalized.estimator(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    # dataa=try(data.generator.and.unpenalized.estimator(f,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv),silent=TRUE)
    # if ('try-error' %in% class(dataa)){
    #   cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
    #   next
    # }else{
    #   dataa=dataa
    # }
    # oracle.fit=FreqID.LT(dataa$Y.oracle, dataa$lin.pred.oracle, dataa$data.oracle, model = "semi-Markov", frailty=frailty, startVals=NULL)
    
    optim_starting_values[which(f==seeds),]=dataa$fit.bp$startingValues[1:11]
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
    
    a=try(GCV.finder.BAR.bp(tol,lambdavec.bar,xi,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,
                            dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,dataa$Y$delta2,dataa$Y$L,
                            dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,
                            dataa$lin.pred,dataa$data, model = "semi-Markov",frailty,dataa$startVals,
                            dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    if ('try-error' %in% class(a)){
      cat("I see a problem in running GCV.finder function for BAR, so, I am jumping to the next iteration \n")
      next
    }else{
      a=a
    }
    
    Y=dataa$Y
    y1=dataa$y1
    y2=dataa$y2
    delta1=dataa$delta1
    delta2=dataa$delta2
    l=dataa$l
    Xmat1=dataa$Xmat1
    Xmat2=dataa$Xmat2
    Xmat3=dataa$Xmat3
    num.of.beta.in.each.risk=length(a$beta.GCV)/3
    b.1=dataa$b.1
    b.2=dataa$b.2
    b.3.y2my1=dataa$b.3.y2my1
    
    #I want to update theta here:
    #based on BAR:
    initial.theta.to.update.bar=a$theta.final
    
    
    logLike.to.update.theta <- function(p) logLike.SCR.SM.LT.bSpline.bp.dropPrevCases.beta.fixed(p,beta1.from.varsel=a$beta.GCV[1:num.of.beta.in.each.risk],
                                                                                                 beta2.from.varsel=a$beta.GCV[(num.of.beta.in.each.risk+1):(2*num.of.beta.in.each.risk)],
                                                                                                 beta3.from.varsel=a$beta.GCV[((2*num.of.beta.in.each.risk)+1):(3*num.of.beta.in.each.risk)] 
                                                                                                 , y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,b.1,  
                                                                                                 b.2,  
                                                                                                 b.3.y2my1)
    optim.control = list(REPORT = 50)
    startVals.bar=initial.theta.to.update.bar#gertting the initial value to optimize from the theta that was estimated in unpenalized method
    fit1 <- optim(startVals.bar, #* runif(length(startVals), 0.9, 1.1), 
                  logLike.to.update.theta, hessian = TRUE, method="Nelder-Mead", control = optim.control)
    BAR.value <- list(estimate=fit1$par, H=fit1$hessian, logLike=-fit1$value, code=fit1$convergence)#, Xmat=list(Xmat1, Xmat2, Xmat3))
    BAR.updated.theta[,which(f==seeds)]=BAR.value$estimate
    
    
    
    
    
    blasso=try(GCV.finder.lasso.bp(lambdavec.lasso,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,
                                   dataa$Y$delta2,dataa$Y$L,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov",
                                   frailty, dataa$startVals,
                                   dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    if ('try-error' %in% class(blasso)){
      cat("I see a problem in running GCV.finder function for lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      blasso=blasso
    }
    
    
    #based on LASSO:
    initial.theta.to.update.lasso=blasso$theta.final
    logLike.to.update.theta <- function(p) logLike.SCR.SM.LT.bSpline.bp.dropPrevCases.beta.fixed(p,beta1.from.varsel=blasso$beta.GCV[1:num.of.beta.in.each.risk],
                                                                                                 beta2.from.varsel=blasso$beta.GCV[(num.of.beta.in.each.risk+1):(2*num.of.beta.in.each.risk)],
                                                                                                 beta3.from.varsel=blasso$beta.GCV[((2*num.of.beta.in.each.risk)+1):(3*num.of.beta.in.each.risk)] 
                                                                                                 , y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,b.1,  
                                                                                                 b.2,  
                                                                                                 b.3.y2my1)
    
    optim.control = list(REPORT = 50)
    startVals.lasso=initial.theta.to.update.lasso#gertting the initial value to optimize from the theta that was estimated in unpenalized method
    fit1 <- optim(startVals.lasso, #* runif(length(startVals), 0.9, 1.1), 
                  logLike.to.update.theta, hessian = TRUE, method="Nelder-Mead", control = optim.control)
    lasso.value <- list(estimate=fit1$par, H=fit1$hessian, logLike=-fit1$value, code=fit1$convergence)#, Xmat=list(Xmat1, Xmat2, Xmat3))
    lasso.updated.theta[,which(f==seeds)]=lasso.value$estimate
    
    
    cada=try(GCV.finder.AdaLasso.bp(lambdavec.ada,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,
                                    dataa$Y$delta2,dataa$Y$L,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov",
                                    frailty, dataa$startVals,
                                    dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    
    # a=try(GCV.finder.lasso(lambdavec,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(cada)){
      cat("I see a problem in running GCV.finder function for adaptive lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      cada=cada
    }
    
    #based on Adaptive LASSO:
    initial.theta.to.update.adalasso=cada$theta.final
    logLike.to.update.theta <- function(p) logLike.SCR.SM.LT.bSpline.bp.dropPrevCases.beta.fixed(p,beta1.from.varsel=cada$beta.GCV[1:num.of.beta.in.each.risk],
                                                                                                 beta2.from.varsel=cada$beta.GCV[(num.of.beta.in.each.risk+1):(2*num.of.beta.in.each.risk)],
                                                                                                 beta3.from.varsel=cada$beta.GCV[((2*num.of.beta.in.each.risk)+1):(3*num.of.beta.in.each.risk)] 
                                                                                                 , y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,b.1,  
                                                                                                 b.2,  
                                                                                                 b.3.y2my1)
    optim.control = list(REPORT = 50)
    startVals.adalasso=initial.theta.to.update.adalasso#gertting the initial value to optimize from the theta that was estimated in unpenalized method
    fit1 <- optim(startVals.adalasso, #* runif(length(startVals), 0.9, 1.1), 
                  logLike.to.update.theta, hessian = TRUE, method="Nelder-Mead", control = optim.control)
    adalasso.value <- list(estimate=fit1$par, H=fit1$hessian, logLike=-fit1$value, code=fit1$convergence)#, Xmat=list(Xmat1, Xmat2, Xmat3))
    adalasso.updated.theta[,which(f==seeds)]=adalasso.value$estimate
    
    # GCV.selected[(which(seed==seeds))]=a$lam.final
    GCV.selected.bar[,(which(f==seeds))]=a$GCV.bar
    selected.beta.based.on.GCV.bar=a$beta.GCV
    TP.bar[(which(f==seeds))]=a$TP.final
    FP.bar[(which(f==seeds))]=a$FP.final
    MSE.bar[(which(f==seeds))]=a$MSE.final
    beta.selected.bar[,(which(f==seeds))]=a$beta.GCV
    FP.bar[(which(f==seeds))]=a$FP.final
    FN.bar[(which(f==seeds))]=a$FN.final
    theta.bar[,(which(f==seeds))]=a$theta.final
    
    GCV.selected.lasso[,(which(f==seeds))]=blasso$GCV.lasso
    selected.beta.based.on.GCV.lasso=blasso$beta.GCV
    TP.lasso[(which(f==seeds))]=blasso$TP.final
    beta.selected.lasso[,(which(f==seeds))]=blasso$beta.GCV
    FP.lasso[(which(f==seeds))]=blasso$FP.final
    MSE.lasso[(which(f==seeds))]=blasso$MSE.final
    FN.lasso[(which(f==seeds))]=blasso$FN.final
    theta.lasso[,(which(f==seeds))]=blasso$theta.final
    
    GCV.selected.ada[,(which(f==seeds))]=cada$GCV.AdaLasso
    selected.beta.based.on.GCV.ada=cada$beta.GCV
    TP.ada[(which(f==seeds))]=cada$TP.final
    beta.selected.ada[,(which(f==seeds))]=cada$beta.GCV
    FP.ada[(which(f==seeds))]=cada$FP.final
    MSE.ada[(which(f==seeds))]=cada$MSE.final
    FN.ada[(which(f==seeds))]=cada$FN.final
    theta.ada[,(which(f==seeds))]=cada$theta.final
    
    
  }
  theta.bar.mean=apply(theta.bar,1,mean)
  theta.ada.mean=apply(theta.ada,1,mean)
  theta.lasso.mean=apply(theta.lasso,1,mean)
  
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
  
  return(list(optim_starting_values=optim_starting_values,BAR.updated.theta=BAR.updated.theta,lasso.updated.theta=lasso.updated.theta,adalasso.updated.theta=adalasso.updated.theta,theta.bar.mean=theta.bar.mean,theta.ada.mean=theta.ada.mean,theta.lasso.mean=theta.lasso.mean,
              theta.ada=theta.ada,theta.lasso=theta.lasso,theta.bar=theta.bar,censoring.rate.final=censoring.rate.final,censoring.observ.rates=censoring.observ.rates,dimension=dimension,dimension.final.mean=dimension.final.mean
              ,SD.bar=SD.bar,beta.selected.bar=beta.selected.bar,GCV.selected.bar=GCV.selected.bar,TP.bar=TP.bar,FP.bar=FP.bar,FN.bar=FN.bar,MC.final.bar=MC.final.bar,MSE.bar=MSE.bar,MSE.calculator.oracle=MSE.calculator.oracle
              ,SD.lasso=SD.lasso,beta.selected.lasso=beta.selected.lasso,GCV.selected.lasso=GCV.selected.lasso,TP.lasso=TP.lasso,FP.lasso=FP.lasso,FN.lasso=FN.lasso,MC.lasso=MC.lasso,MC.final.lasso=MC.final.lasso,MSE.lasso=MSE.lasso
              ,SD.ada=SD.ada,beta.selected.ada=beta.selected.ada,GCV.selected.ada=GCV.selected.ada,TP.ada=TP.ada,FP.ada=FP.ada,FN.ada=FN.ada,MC.ada=MC.ada,MC.final.ada=MC.final.ada,MSE.ada=MSE.ada,SD.oracle=SD.oracle))
}



repeat.get.last.result.bp=function(seeds,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,lambdavec.bar,lambdavec.lasso,lambdavec.ada,theta.given,frailty_distribution){
  for (f in seeds){
    cat("this is risks.stat:",risks.status,"\n")
    cat("Now running seed",f,"\n")
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method)
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status)
    if (length(beta1.true)==12){
      dataa=try(data.generator.and.unpenalized.estimator.12cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if (length(beta1.true)==15){
      dataa=try(data.generator.and.unpenalized.estimator.15cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if (length(beta1.true)==16){
      dataa=try(data.generator.and.unpenalized.estimator.16cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    
    if (length(beta1.true)==20){
      dataa=try(data.generator.and.unpenalized.estimator.20cov(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    if (length(beta1.true)==10){
      dataa=try(data.generator.and.unpenalized.estimator(f,mu1,mu2,sigma1,sigma2,p1,p2,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv,theta.given,frailty_distribution),silent=TRUE)
      if ('try-error' %in% class(dataa)){
        cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
        next
      }else{
        dataa=dataa
      }
    }
    # dataa=try(data.generator.and.unpenalized.estimator(f,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv),silent=TRUE)
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
    
    # oracle.fit=try(FreqID.LT(dataa$Y.oracle.1, dataa$lin.pred.oracle, dataa$data.oracle.1, model = "semi-Markov", frailty=frailty, startVals=NULL,method),silent=TRUE)
    # if ('try-error' %in% class(oracle.fit)){
    #   cat("I see a problem in running FreqID.LT function, so, I am jumping to the next iteration \n")
    #   next
    # }else{
    #   oracle.fit=oracle.fit
    # }
    MSE.calculator.oracle[which(f==seeds)]=t(as.matrix((dataa$unpen.est.oracle)-dataa$beta.truth.oracle))%*%(dataa$E_ZZprime.oracle)%*%as.matrix((dataa$unpen.est.oracle)-dataa$beta.truth.oracle)
    # a=GCV.finder.BAR(tol,lambdavec,xi,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL)
    
    a=try(GCV.finder.BAR.bp(tol,lambdavec.bar,xi,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,
                            dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,dataa$Y$delta2,dataa$Y$L,
                            dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,
                            dataa$lin.pred,dataa$data, model = "semi-Markov",frailty,dataa$startVals,
                            dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    if ('try-error' %in% class(a)){
      cat("I see a problem in running GCV.finder function for BAR, so, I am jumping to the next iteration \n")
      next
    }else{
      a=a
    }
    
    blasso=try(GCV.finder.lasso.bp(lambdavec.lasso,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,
                                   dataa$Y$delta2,dataa$Y$L,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov",
                                   frailty, dataa$startVals,
                                   dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    if ('try-error' %in% class(blasso)){
      cat("I see a problem in running GCV.finder function for lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      blasso=blasso
    }
    
    cada=try(GCV.finder.AdaLasso.bp(lambdavec.ada,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,
                                    dataa$Y$delta2,dataa$Y$L,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov",
                                    frailty, dataa$startVals,
                                    dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    
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
    theta.bar[,(which(f==seeds))]=a$theta.final
    
    GCV.selected.lasso[,(which(f==seeds))]=blasso$GCV.lasso
    selected.beta.based.on.GCV.lasso=blasso$beta.GCV
    TP.lasso[(which(f==seeds))]=blasso$TP.final
    beta.selected.lasso[,(which(f==seeds))]=blasso$beta.GCV
    FP.lasso[(which(f==seeds))]=blasso$FP.final
    MSE.lasso[(which(f==seeds))]=blasso$MSE.final
    FN.lasso[(which(f==seeds))]=blasso$FN.final
    theta.lasso[,(which(f==seeds))]=blasso$theta.final
    
    GCV.selected.ada[,(which(f==seeds))]=cada$GCV.AdaLasso
    selected.beta.based.on.GCV.ada=cada$beta.GCV
    TP.ada[(which(f==seeds))]=cada$TP.final
    beta.selected.ada[,(which(f==seeds))]=cada$beta.GCV
    FP.ada[(which(f==seeds))]=cada$FP.final
    MSE.ada[(which(f==seeds))]=cada$MSE.final
    FN.ada[(which(f==seeds))]=cada$FN.final
    theta.ada[,(which(f==seeds))]=cada$theta.final
    
    
  }
  theta.bar.mean=apply(theta.bar,1,mean)
  theta.ada.mean=apply(theta.ada,1,mean)
  theta.lasso.mean=apply(theta.lasso,1,mean)
  
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
  
  return(list(theta.bar.mean=theta.bar.mean,theta.ada.mean=theta.ada.mean,theta.lasso.mean=theta.lasso.mean,
              theta.ada=theta.ada,theta.lasso=theta.lasso,theta.bar=theta.bar,censoring.rate.final=censoring.rate.final,censoring.observ.rates=censoring.observ.rates,dimension=dimension,dimension.final.mean=dimension.final.mean
              ,SD.bar=SD.bar,beta.selected.bar=beta.selected.bar,GCV.selected.bar=GCV.selected.bar,TP.bar=TP.bar,FP.bar=FP.bar,FN.bar=FN.bar,MC.final.bar=MC.final.bar,MSE.bar=MSE.bar,MSE.calculator.oracle=MSE.calculator.oracle
              ,SD.lasso=SD.lasso,beta.selected.lasso=beta.selected.lasso,GCV.selected.lasso=GCV.selected.lasso,TP.lasso=TP.lasso,FP.lasso=FP.lasso,FN.lasso=FN.lasso,MC.lasso=MC.lasso,MC.final.lasso=MC.final.lasso,MSE.lasso=MSE.lasso
              ,SD.ada=SD.ada,beta.selected.ada=beta.selected.ada,GCV.selected.ada=GCV.selected.ada,TP.ada=TP.ada,FP.ada=FP.ada,FN.ada=FN.ada,MC.ada=MC.ada,MC.final.ada=MC.final.ada,MSE.ada=MSE.ada,SD.oracle=SD.oracle))
}

repeat.get.last.result.grp=function(seeds,tol,n,m,yvalslength,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,weibull.param.log,lambdavec.lasso,lambdavec.ada,lambdavec.bar,xi,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli){
  for (f in seeds){
    cat("this is risks.stat:",risks.status,"\n")
    cat("Now running seed",f,"\n")
    # dataa=data.generator(f,n,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1=G1,G2=G2,G3=G3,G4=G4,G1.distr=G1.distr,G2.distr=G2.distr,G3.distr=G3.distr,G4.distr=G4.distr,mu.bernoulli=mu.bernoulli,grp.or.indiv)
    dataa=try(data.generator.and.unpenalized.estimator(f,n,m,yvalslength,weibull.param.log,c1,c2,lt1,lt2,beta1.true,beta2.true,beta3.true,rho,method,risks.status,G1,G2,G3,G4,G1.distr,G2.distr,G3.distr,G4.distr,mu.bernoulli,grp.or.indiv),silent=TRUE)
    if ('try-error' %in% class(dataa)){
      cat("I see a problem in running data.generetor function, so, I am jumping to the next iteration \n")
      next
    }else{
      dataa=dataa
    }
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
    
    a=try(GCV.finder.BAR.bp.grp(G1,G2,G3,G4,tol,lambdavec.bar,xi,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,
                                dataa$Y,dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,dataa$Y$delta2,dataa$Y$L,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,
                                dataa$lin.pred,dataa$data, model = "semi-Markov",frailty, dataa$startVals,dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    
    
    if ('try-error' %in% class(a)){
      cat("I see a problem in running GCV.finder function, so, I am jumping to the next iteration \n")
      next
    }else{
      a=a
    }
    
    # blasso=GCV.finder.lasso.grp(G1,G2,G3,G4,lambdavec.lasso,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL)
    
    blasso=try(GCV.finder.lasso.bp.grp(G1,G2,G3,G4,lambdavec.lasso,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,
                                       dataa$Y,dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,dataa$Y$delta2,dataa$Y$L,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,
                                       dataa$lin.pred,dataa$data, model = "semi-Markov",frailty, dataa$startVals,dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    # blasso=try(GCV.finder.lasso.grp(G1,G2,G3,G4,lambdavec.lasso,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL),silent=TRUE)
    if ('try-error' %in% class(blasso)){
      cat("I see a problem in running GCV.finder function for lasso, so, I am jumping to the next iteration \n")
      next
    }else{
      blasso=blasso
    }
    
    # cada=GCV.finder.AdaLasso.grp(G1,G2,G3,G4,lambdavec.ada,dataa$weibull.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,dataa$Y,dataa$y1,dataa$y2,dataa$delta1,dataa$delta2,dataa$l,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,dataa$lin.pred,dataa$data, model = "semi-Markov", frailty=TRUE, startVals=NULL)
    cada=try(GCV.finder.AdaLasso.bp.grp(G1,G2,G3,G4,lambdavec.lasso,dataa$nuisance.estimates,dataa$beta.truth,dataa$E_ZZprime,dataa$unpen.est,
                                        dataa$Y,dataa$Y$y1,dataa$Y$y2,dataa$Y$delta1,dataa$Y$delta2,dataa$Y$L,dataa$Xmat1,dataa$Xmat2,dataa$Xmat3,
                                        dataa$lin.pred,dataa$data, model = "semi-Markov",frailty, dataa$startVals,dataa$b.1,dataa$b.2,dataa$b.3.y2my1),silent=TRUE)
    
    
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
  G.total.bar=(0.2*G1.percentage.f.bar)+(0.2*G2.percentage.f.bar)+(0.3*G3.percentage.f.bar)+(0.3*G4.percentage.f.bar)
  
  
  G1.percentage.f.lasso=sum(G1.percentage.lasso[!is.na(G1.percentage.lasso)])/length(G1.percentage.lasso[!is.na(G1.percentage.lasso)])
  G2.percentage.f.lasso=sum(G2.percentage.lasso[!is.na(G2.percentage.lasso)])/length(G2.percentage.lasso[!is.na(G2.percentage.lasso)])
  G3.percentage.f.lasso=sum(G3.percentage.lasso[!is.na(G3.percentage.lasso)])/length(G3.percentage.lasso[!is.na(G3.percentage.lasso)])
  G4.percentage.f.lasso=sum(G4.percentage.lasso[!is.na(G4.percentage.lasso)])/length(G4.percentage.lasso[!is.na(G4.percentage.lasso)])
  G.total.lasso=(0.2*G1.percentage.f.lasso)+(0.2*G2.percentage.f.lasso)+(0.3*G3.percentage.f.lasso)+(0.3*G4.percentage.f.lasso)
  
  
  G1.percentage.f.ada=sum(G1.percentage.ada[!is.na(G1.percentage.ada)])/length(G1.percentage.ada[!is.na(G1.percentage.ada)])
  G2.percentage.f.ada=sum(G2.percentage.ada[!is.na(G2.percentage.ada)])/length(G2.percentage.ada[!is.na(G2.percentage.ada)])
  G3.percentage.f.ada=sum(G3.percentage.ada[!is.na(G3.percentage.ada)])/length(G3.percentage.ada[!is.na(G3.percentage.ada)])
  G4.percentage.f.ada=sum(G4.percentage.ada[!is.na(G4.percentage.ada)])/length(G4.percentage.ada[!is.na(G4.percentage.ada)])
  G.total.ada=(0.2*G1.percentage.f.ada)+(0.2*G2.percentage.f.ada)+(0.3*G3.percentage.f.ada)+(0.3*G4.percentage.f.ada)
  
  return(list(G.total.bar=G.total.bar,G1.percentage.f.bar=G1.percentage.f.bar,G2.percentage.f.bar=G2.percentage.f.bar,G3.percentage.f.bar=G3.percentage.f.bar,G4.percentage.f.bar=G4.percentage.f.bar,
              G.total.lasso=G.total.lasso,G1.percentage.f.lasso=G1.percentage.f.lasso,G2.percentage.f.lasso=G2.percentage.f.lasso,G3.percentage.f.lasso=G3.percentage.f.lasso,G4.percentage.f.lasso=G4.percentage.f.lasso,
              G.total.ada=G.total.ada,G1.percentage.f.ada=G1.percentage.f.ada,G2.percentage.f.ada=G2.percentage.f.ada,G3.percentage.f.ada=G3.percentage.f.ada,G4.percentage.f.ada=G4.percentage.f.ada,
              censoring.rate.final=censoring.rate.final,censoring.observ.rates=censoring.observ.rates,dimension=dimension
              ,SD.bar=SD.bar,beta.selected.bar=beta.selected.bar,GCV.selected.bar=GCV.selected.bar,TP.bar=TP.bar,FP.bar=FP.bar,FN.bar=FN.bar,MC.final.bar=MC.final.bar,MSE.bar=MSE.bar,MSE.calculator.oracle=MSE.calculator.oracle
              ,SD.lasso=SD.lasso,beta.selected.lasso=beta.selected.lasso,GCV.selected.lasso=GCV.selected.lasso,TP.lasso=TP.lasso,FP.lasso=FP.lasso,FN.lasso=FN.lasso,MC.lasso=MC.lasso,MC.final.lasso=MC.final.lasso,MSE.lasso=MSE.lasso
              ,SD.ada=SD.ada,beta.selected.ada=beta.selected.ada,GCV.selected.ada=GCV.selected.ada,TP.ada=TP.ada,FP.ada=FP.ada,FN.ada=FN.ada,MC.ada=MC.ada,MC.final.ada=MC.final.ada,MSE.ada=MSE.ada,SD.oracle=SD.oracle
  ))
}












