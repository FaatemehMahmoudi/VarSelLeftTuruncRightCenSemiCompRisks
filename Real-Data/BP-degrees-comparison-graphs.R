m1=2
m2=2
m3=3
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

a=c(-8.839995308 , -9.152434611 ,-15.707705064)
a1=c(-9.520910683 , -8.650990264, -10.429091399 ,-10.534416876,-13.181870250 ,-7.226307337,-25.289015778)
b=c(-12.270242526 ,-10.796645677 ,-15.453769578 )
b1=c(-11.701442737, -13.727111643 ,-11.967549923 , -9.443884958 ,-13.487411428 ,-7.343357993,-19.593029616)
c= c(-8.394288258 , -7.336739345  ,-8.391945016 ,-15.163452368)
c1=c(-8.744380222 ,-7.800910824 ,-8.365696790 ,-4.728466706 ,-13.698021779,-8.120005639,-20.805236712)
t1=seq(min(y1),max(y1),0.1)
t2=seq(min(y2),max(y2),0.1)

length(t1)



#first trans:
plot.y1.bp.m223=exp(as.vector(matrix(fit.bh.crc.regm223$estimate[1:3], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))
plot.y1.bp.m334=exp(as.vector(matrix(fit.bh.crc.regm334$estimate[1:4], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))
plot.y1.bp.m666=exp(as.vector(matrix(fit.bh.crc.regm666$estimate[1:7], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))
plot.y1.bp.m999=exp(as.vector(matrix(fit.bh.crc.regm999$estimate[1:10], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))
#baseline haz:
plot(t1,plot.y1.bp.m223,type="l",col="red")
lines(t1,plot.y1.bp.m334)
lines(t1,plot.y1.bp.m666,col="blue")
lines(t1,plot.y1.bp.m999,col="green")
#cum baseline haz:
plot(t1,cumsum(plot.y1.bp.m223),type="l",col="red")
lines(t1,cumsum(plot.y1.bp.m334))
lines(t1,cumsum(plot.y1.bp.m666),col="blue")
lines(t1,cumsum(plot.y1.bp.m999),col="green")
legend(800, 0.25, legend=c("m=(2,2,3)", "m=(3,3,4)","m=c(6,6,6)","m=c(9,9,9)"),
       col=c("red", "black","blue","green"),lty=c(1,1,1,1),cex=0.3)


#2ns trans:
plot.y2.bp.m223=exp(as.vector(matrix(fit.bh.crc.regm223$estimate[4:6], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
plot.y2.bp.m666=exp(as.vector(matrix(fit.bh.crc.regm666$estimate[8:14], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
# plot.y2.bp.m667=exp(as.vector(matrix(fit.bh.crc.regm667$estimate[8:14], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
plot.y2.bp.m334=exp(as.vector(matrix(fit.bh.crc.regm334$estimate[5:8], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
plot.y2.bp.m999=exp(as.vector(matrix(fit.bh.crc.regm999$estimate[11:20], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
#baseline haz:
plot(t2,plot.y2.bp.m223,type="l",col="red",ylim=c(0,40e-06))
lines(t2,plot.y2.bp.m666,col="blue")
lines(t2,plot.y2.bp.m334)
lines(t2,plot.y2.bp.m999,col="green")
#cum baseline haz:
plot(t2,cumsum(plot.y2.bp.m223),col="red",type="l",ylim=c(0,0.7))
lines(t2,cumsum(plot.y2.bp.m666),type="l",col="blue")
lines(t2,cumsum(plot.y2.bp.m334))
lines(t2,cumsum(plot.y2.bp.m999),col="green")
legend(100, 0.6, legend=c("m=(2,2,3)", "m=(3,3,4)","m=c(6,6,6)","m=c(9,9,9)"),
       col=c("red", "black","blue","green"),lty=c(1,1,1,1),cex=0.3)


#3rd trans:
t3=seq(min(y2-y1),max(y2-y1),0.1)
plot.y3.bp.m223=exp(as.vector(matrix(fit.bh.crc.regm223$estimate[7:10] , nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
plot.y3.bp.m334=exp(as.vector(matrix(fit.bh.crc.regm334$estimate[9:13], nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
plot.y3.bp.m666=exp(as.vector(matrix(fit.bh.crc.regm666$estimate[15:21], nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
plot.y3.bp.m999=exp(as.vector(matrix(fit.bh.crc.regm999$estimate[21:30], nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
#baseline haz:
plot(t3,plot.y3.bp.m223,type="l",col="red",ylim=c(0,0.00039))
lines(t3,plot.y3.bp.m334)
lines(t3,plot.y3.bp.m666,col="blue")
lines(t3,plot.y3.bp.m999,col="green")
#cum baseline haz:
plot(t3,cumsum(plot.y3.bp.m223),type="l",col="red",ylim=c(0,4.7))
lines(t3,cumsum(plot.y3.bp.m334))
lines(t3,cumsum(plot.y3.bp.m666),col="blue")
lines(t3,cumsum(plot.y3.bp.m999),col="green")
legend(1000, 1, legend=c("m=(2,2,3)", "m=(3,3,4)","m=c(6,6,6)","m=c(9,9,9)"),
       col=c("red", "black","blue","green"),lty=c(1,1,1,1),cex=0.3)

filename=c("m223-m334-m666-m999")
sink(filename)
print("m223:")
print(fit.bh.crc.regm223)

print("\n\n\n")
print("m334:")
print(fit.bh.crc.regm334)
print("\n\n\n")
print("m666:")
print(fit.bh.crc.regm666)
print("\n\n\n")
print("m999:")
print(fit.bh.crc.regm999)
sink()









#######################################################
BIC.calc=function(m,param){
  m1=m[1]
  m2=m[2]
  m3=m[3]
  degree.param=length(param)
  
  y1=dat.crc.reg$y1
  y2=dat.crc.reg$y2
  delta1=dat.crc.reg$delta1
  delta2=dat.crc.reg$delta2
  lin.pred=lin.pred.crc
  data=dat.crc.reg
  Xmat1=as.matrix(data[,(6:17)])
  Xmat2=Xmat1
  Xmat3=Xmat1
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
  
  #second trans:
  b.2.event.bp=bernsteinPoly(x=y2[delta2==1],degree=m2,intercept = T,Boundary.knots = bdy.knots.b.2)
  b.2.bp=predict(b.2.event.bp,y2)
  
  #third trans:
  b.3.event.bp=bernsteinPoly(x=(y2-y1)[delta1==1],degree = m3,intercept = T,Boundary.knots = bdy.knots.b.3.y2my1)
  b.3.y2my1.bp=predict(b.3.event.bp,y2-y1)
  loglik=logLike.SCR.SM.LT.bSpline.bp.dropPrevCases(param, y1, y2, delta1, delta2, data$L, Xmat1, Xmat2, Xmat3, frailty=TRUE,
                                                                b.1.bp,  
                                                                b.2.bp,  
                                                                b.3.y2my1.bp)
  BIC=-2*loglik
  
}



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




#AIC, BIC:
#m=c(2,2,3):
AICm223=(2*fit.bh.crc.regm223$logLike)+(2*length(fit.bh.crc.regm223$estimate))
AICm334=(2*fit.bh.crc.regm334$logLike)+(2*length(fit.bh.crc.regm334$estimate))
AICm666=(2*fit.bh.crc.regm666$logLike)+(2*length(fit.bh.crc.regm666$estimate))
AICm999=(2*fit.bh.crc.regm999$logLike)+(2*length(fit.bh.crc.regm999$estimate))


#calculaytion based on what DR. Lu mentioned:
BICm223=((2*fit.bh.crc.regm223$logLike)+(log(n)*length(fit.bh.crc.regm223$estimate)))
BICm334=((2*fit.bh.crc.regm334$logLike)+(log(n)*length(fit.bh.crc.regm334$estimate)))
BICm666=((2*fit.bh.crc.regm666$logLike)+(log(n)*length(fit.bh.crc.regm666$estimate)))
BICm999=((2*fit.bh.crc.regm999$logLike)+(log(n)*length(fit.bh.crc.regm999$estimate)))

min(BICm223,BICm334,BICm666,BICm999)


sink("filename")
cat("AIC's for 223,334,666,999:",AICm223,"\n",AICm334,"\n",AICm666,"\n",AICm999,"\n")
cat("BIC's for 223,334,666,999:",BICm223,"\n",BICm334,"\n",BICm666,"\n",BICm999,"\n")
sink()





##############################
#m=c(2,2,3):
est.m223=scan()
# -8.839995308  -9.152434611 -15.707705064 -12.270242526 -10.796645677 -15.453769578 
# -8.394288258  -7.336739345  -8.391945016 -15.163452368  -2.402022783   0.104066947 
# -0.391087615  -0.080279248  -0.009020569   0.261586258   0.229130312   0.252200334 
# 0.050786534   0.260034968   0.396359907   0.349414335   0.563771745   0.226492830 
# 0.152439292   0.035328786   0.009071025   0.356905397  -0.711071749   0.178874276 
# -0.032188130   0.372827483   0.026699255   0.423358983   0.871770855   0.176626892 
# 0.330629965   0.300942396   0.009573824   0.265657864  -0.575287797   0.169492763 
# 0.020304043   0.042546235   0.183281879   0.042564770   0.520873185 
est.m334=scan()
# -9.212414206  -8.954124973 -12.179982188 -16.278724207 -12.206570556 -11.769896156 
# -8.509985007 -19.544486755  -8.882228552  -7.842991015  -8.051881353 -11.423948165 
# -14.840841678  -1.194156742  -0.151815191  -0.590568171  -0.129019245  -0.004530600 
# 0.281795404   0.110203000   0.339831441   0.043593611   0.229631986   0.529991405 
# 0.215459150   0.594747947  -0.114554029  -0.175062799  -0.173816587   0.003932156 
# 0.669531215   0.189589388   0.503001579  -0.138904348   0.394557667   0.199851927 
# 0.294643725   1.609214181   0.128551879   0.337299683   0.199782101   0.012139385 
# 0.358222703  -0.457446726   0.265800750   0.040571014   0.044142773   0.249726674 
# 0.122542188   0.416436623 
est.m666=scan()
# -9.520910683  -8.650990264 -10.429091399 -10.534416876 -13.181870250  -7.226307337 
# -25.289015778 -11.701442737 -13.727111643 -11.967549923  -9.443884958 -13.487411428 
# -7.343357993 -19.593029616  -8.744380222  -7.800910824  -8.365696790  -4.728466706 
# -13.698021779  -8.120005639 -20.805236712  -0.838185614  -0.166042445  -0.547028911 
# -0.193144320  -0.001445752   0.303197073   0.289985764   0.192562448   0.053413086 
# 0.233834699   0.524359453   0.260037812   0.562308674  -0.226784521   0.023922453 
# 0.038383520   0.002105782   1.065352488  -1.821427547   0.729681726  -0.086634715 
# 0.030699609   0.358526594   0.501406531   1.339201781   0.101097342   0.359207019 
# 0.020441761   0.010476885   0.447148922  -0.240049914   0.187848609   0.034831165 
# 0.073815634   0.216235218   0.025423406   0.430838319 
est.m999=scan()
# -9.289623251  -8.520482558  -9.127053625 -10.805569239  -6.176722912 -24.306637934
# 15.220014588 -50.766161872  35.532760140 -73.836672962 -11.079241665 -10.809211402
# -11.118101319 -12.457468859  -6.248718667 -18.976565766   7.380839098 -33.442462221
# 15.382351005 -59.801173903  -8.605017741  -8.150473822  -7.083215045  -9.140888348
# -3.541201113 -18.748713720   5.055835266 -33.684243233  17.583178963 -44.873466331
# -1.468442045   0.001982091  -0.483901238  -0.135048393  -0.004227175   0.284645208
# 0.170856522   0.261173404   0.044011392   0.188160383   0.453720855   0.302115936
# 0.595127799  -0.364140967  -0.361505868   0.021599218  -0.009681376   0.170811595
# -0.644153110   0.057369466  -0.019871905   0.013220309   0.378172702   0.298685658
# 0.684628621   0.229111425   0.412102737   0.202637277   0.008461245   0.376987550
# -0.315916687   0.244651020   0.024494150   0.066175550   0.182401147   0.113064043
# 0.523422797









t1=seq(0,max(y1),0.1)
t2=t1
t3=seq(0,max(y2-y1),0.1)
#first trans:
plot.y1.bp.m223.haz=exp(as.vector(matrix(est.m223[1:3], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))*exp(t(est.m223[(3+3+4+1+1):length(est.m223)])%*%t(as.matrix(cbind(dat.crc.reg[,(6:17)],dat.crc.reg[,(6:17)],dat.crc.reg[,(6:17)]))))
plot.y1.bp.m334.haz=exp(as.vector(matrix(fit.bh.crc.regm334$estimate[1:4], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))
plot.y1.bp.m666.haz=exp(as.vector(matrix(fit.bh.crc.regm666$estimate[1:7], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))
plot.y1.bp.m999.haz=exp(as.vector(matrix(fit.bh.crc.regm999$estimate[1:10], nrow=1) %*% t(as.matrix(predict(b.1.bp , t1)))))
#baseline haz:
plot(t1,plot.y1.bp.m223,type="l",col="red")
lines(t1,plot.y1.bp.m334)
lines(t1,plot.y1.bp.m666,col="blue")
lines(t1,plot.y1.bp.m999,col="green")
#cum baseline haz:
plot(t1,cumsum(plot.y1.bp.m223),type="l",col="red")
lines(t1,cumsum(plot.y1.bp.m334))
lines(t1,cumsum(plot.y1.bp.m666),col="blue")
lines(t1,cumsum(plot.y1.bp.m999),col="green")
legend(800, 0.25, legend=c("m=(2,2,3)", "m=(3,3,4)","m=c(6,6,6)","m=c(9,9,9)"),
       col=c("red", "black","blue","green"),lty=c(1,1,1,1),cex=0.3)


#2ns trans:
plot.y2.bp.m223=exp(as.vector(matrix(fit.bh.crc.regm223$estimate[4:6], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
plot.y2.bp.m666=exp(as.vector(matrix(fit.bh.crc.regm666$estimate[8:14], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
# plot.y2.bp.m667=exp(as.vector(matrix(fit.bh.crc.regm667$estimate[8:14], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
plot.y2.bp.m334=exp(as.vector(matrix(fit.bh.crc.regm334$estimate[5:8], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
plot.y2.bp.m999=exp(as.vector(matrix(fit.bh.crc.regm999$estimate[11:20], nrow=1) %*% t(as.matrix(predict(b.2.bp , t2)))))
#baseline haz:
plot(t2,plot.y2.bp.m223,type="l",col="red",ylim=c(0,40e-06))
lines(t2,plot.y2.bp.m666,col="blue")
lines(t2,plot.y2.bp.m334)
lines(t2,plot.y2.bp.m999,col="green")
#cum baseline haz:
plot(t2,cumsum(plot.y2.bp.m223),col="red",type="l",ylim=c(0,0.7))
lines(t2,cumsum(plot.y2.bp.m666),type="l",col="blue")
lines(t2,cumsum(plot.y2.bp.m334))
lines(t2,cumsum(plot.y2.bp.m999),col="green")
legend(100, 0.6, legend=c("m=(2,2,3)", "m=(3,3,4)","m=c(6,6,6)","m=c(9,9,9)"),
       col=c("red", "black","blue","green"),lty=c(1,1,1,1),cex=0.3)


#3rd trans:
t3=seq(min(y2-y1),max(y2-y1),0.1)
plot.y3.bp.m223=exp(as.vector(matrix(fit.bh.crc.regm223$estimate[7:10] , nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
plot.y3.bp.m334=exp(as.vector(matrix(fit.bh.crc.regm334$estimate[9:13], nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
plot.y3.bp.m666=exp(as.vector(matrix(fit.bh.crc.regm666$estimate[15:21], nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
plot.y3.bp.m999=exp(as.vector(matrix(fit.bh.crc.regm999$estimate[21:30], nrow=1) %*% t(as.matrix(predict(b.3.y2my1.bp , t3)))))
#baseline haz:
plot(t3,plot.y3.bp.m223,type="l",col="red",ylim=c(0,0.00039))
lines(t3,plot.y3.bp.m334)
lines(t3,plot.y3.bp.m666,col="blue")
lines(t3,plot.y3.bp.m999,col="green")
#cum baseline haz:
plot(t3,cumsum(plot.y3.bp.m223),type="l",col="red",ylim=c(0,4.7))
lines(t3,cumsum(plot.y3.bp.m334))
lines(t3,cumsum(plot.y3.bp.m666),col="blue")
lines(t3,cumsum(plot.y3.bp.m999),col="green")
legend(1000, 1, legend=c("m=(2,2,3)", "m=(3,3,4)","m=c(6,6,6)","m=c(9,9,9)"),
       col=c("red", "black","blue","green"),lty=c(1,1,1,1),cex=0.3)
