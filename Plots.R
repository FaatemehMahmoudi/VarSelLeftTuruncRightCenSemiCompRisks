#this file is for creating plots of baseline hazard functions along with 
#the cumulative baseline hazard functions for Weibull, and Bernstein Polynomials
#with two sets of parameters. 

#if you are running this priogram, you need to either read the files already written using the 
#"get-baseline-hazards-points.R" file, or run that programn first to create files yourself. 

#If you read the files, you are ready to use this file to create the plots. 

rm(list=ls())
setwd("/Users/fatemehmahmoudi/Desktop/Desktop/Project2/Codes-bin/Plots/DropBox-R-codes-2nd-project/Plots")
source("functions.R")
beta1.true = c(-0.9,-0.8,0.8,0.00,0.00)
beta2.true = c(-0.9,-0.8,0.8,0.8,0.00)
beta3.true = c(0.9,0.8,0.8,0.9,0.00)

#which set of parameters we are using: Ours or the ones close to Lee et al.'s paper setting?
# param.set="ours"
param.set="original"
method="optim"
frailty=TRUE
if (param.set=="ours"){
  t1 = t2 = seq(0.01,45,0.05); t3 = seq(0.01,7.5,0.05)
  lefttrunc=7.5
}else{
  t1 = t2 = seq(0.01,40,0.05); t3 = seq(0.01,8,0.05)
  lefttrunc=8
}
#rep should be in this format as it represents two things:
#1-number of replications/iterations, and 2-feeds the data generation function with seed.
rep=1:199
#BP degree:
# m=c(6,6,6)
m=c(2,2,3)
m1=m[1]
m2=m[2]
m3=m[3]
yvalslength=500
if (param.set=="ours"){
  c=45
}else{
  c=40
}
if (param.set=="ours"){
  WB.simData <-genWB.simData.ours(n = 500, frailty = T, seed = 1)
}else{
  WB.simData <-genWB.simData(n = 500, frailty = T, seed = 1)
}


# Read BP files -----------------------------------------------------------
#Our params:
#m223:
setwd("/Users/fatemehmahmoudi/Desktop/Desktop/Project2/Codes-bin/Plots/DropBox-R-codes-2nd-project/Plots/BH data/close-to-Lee-params/m666")
y.axis.plot.1.b.h=read.csv("y-axis-trans-1m666-200rep")
y.axis.plot.2.b.h=read.csv("y-axis-trans-2m666-200rep")
y.axis.plot.3.b.h=read.csv("y-axis-trans-3m666-200rep")
#Close to Lee params:
# y.axis.plot.1.b.h=read.csv()
# y.axis.plot.2.b.h=read.csv()
# y.axis.plot.3.b.h=read.csv()
# Read Weibull Files ------------------------------------------------------
#our params:
# setwd("/Users/fatemehmahmoudi/Desktop/Desktop/Project2/Codes-bin/Plots/DropBox-R-codes-2nd-project/Plots/Weibull data/our-params/Baseline-haz")
# y.axis.plot.1.w.h=read.csv("y-axis-trans-1.w.h.txt")
# y.axis.plot.2.w.h=read.csv("y-axis-trans-2.w.h.txt")
# y.axis.plot.3.w.h=read.csv("y-axis-trans-3.w.h.txt")
# #set the path to where the files of Weibull-->our params-->cumulative baseline hazards are stored:
# setwd("/Users/fatemehmahmoudi/Desktop/Desktop/Project2/Codes-bin/Plots/DropBox-R-codes-2nd-project/Plots/Weibull data/our-params/Cum-Baseline-HAZ")
# y.axis.plot.1.w.H=read.csv("y-axis-trans-1.w.cum.txt")
# y.axis.plot.2.w.H=read.csv("y-axis-trans-2.w.cum.txt")
# y.axis.plot.3.w.H=read.csv("y-axis-trans-3.w.cum.txt")
#close to Lee params:
#set the path to where the files of Weibull-->close to Lee et al.'s params-->baseline hazards are stored:
#setwd("/Users/fatemehmahmoudi/Desktop/Desktop/Project2/Codes-bin/Plots/DropBox-R-codes-2nd-project/Plots/Weibull data/close-to-Lee-params/baseline-haz")
# y.axis.plot.1.w.h=read.csv("y-axis-trans-1.w.h.txt")
# y.axis.plot.2.w.h=read.csv("y-axis-trans-1.w.h.txt")
# y.axis.plot.3.w.h=read.csv("y-axis-trans-1.w.h.txt")
#set the path to where the files of Weibull-->close to Lee et al.'s params-->cumulative baseline hazards are stored:
#setwd("/Users/fatemehmahmoudi/Desktop/Desktop/Project2/Codes-bin/Plots/DropBox-R-codes-2nd-project/Plots/Weibull data/close-to-Lee-params/cum-baseline-haz")
# y.axis.plot.1.w.H=read.csv("y-axis-trans-1.w.cum.txt")
# y.axis.plot.2.w.H=read.csv("y-axis-trans-1.w.cum.txt")
# y.axis.plot.3.w.H=read.csv("y-axis-trans-1.w.cum.txt")

#run this section if you want to create plots of Weibull baseline hazards:
# Weibull Plots -----------------------------------------------------------
apply(y.axis.plot.1.w.h,1,mean)->mean.1.h
apply(y.axis.plot.2.w.h,1,mean)->mean.2.h
apply(y.axis.plot.3.w.h,1,mean)->mean.3.h
apply(y.axis.plot.1.w.H,1,mean)->mean.1.H
apply(y.axis.plot.2.w.H,1,mean)->mean.2.H
apply(y.axis.plot.3.w.H,1,mean)->mean.3.H
#1-st trans:
#base haz:
par(mfrow = c(2,1), mar=c(4.1, 4.1, 1.1, 1.1))
plot(WB.simData$plot.WB1, type="l", xlim=c(0,c), ylim=c(0,0.15), ylab=expression(lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
for (i in rep){
  lines(t1,y.axis.plot.1.w.h[,i],col="cyan2",lty="dotted")
}
lines(WB.simData$plot.WB1, type="l", ylab=expression(lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
lines(t1, mean.1.h,col="red2",lwd=2)
legend(0.1, 0.14, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan2","red2"),lty=c(1,3,1),cex=0.3,bty="n")

#cum base haz:
#ours:
#plot(WB.simData$plot.WB1.H$x,WB.simData$plot.WB1.H$y, type="l", xlim=c(0,c), ylim=c(0,4), ylab=expression(Lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
#orig:
plot(WB.simData$plot.WB1.H$x,WB.simData$plot.WB1.H$y, type="l", xlim=c(0,c), ylim=c(0,2.5), ylab=expression(Lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
for (i in rep){
  lines(t1,y.axis.plot.1.w.H[,i] ,col="cyan4",lty="dotted")
}
lines(WB.simData$plot.WB1.H, type="l", ylab=expression(lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
lines(t1, mean.1.H,col="red2",lwd=2)
#ours:
# legend(0.1, 3.8, legend=c("truth", "estimates","mean of estimates"),col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.32,bty="n")
#orig:
legend(0.1, 2.4, legend=c("truth", "estimates","mean of estimates"),col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.32,bty="n")

#2-nd trans:
#base haz:
par(mfrow = c(2,1), mar=c(4.1, 4.1, 1.1, 1.1))
plot(WB.simData$plot.WB2, type="l", xlim=c(0,c), ylim=c(0,0.17), ylab=expression(lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
for (i in rep){
  lines(t2, y.axis.plot.2.w.h[,i],col="cyan2",lty="dotted")
}
lines(WB.simData$plot.WB2, type="l", ylab=expression(lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
lines(t2, mean.2.h,col="red2",lwd=2)
legend(0.1, 0.16, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan2","red2"),lty=c(1,3,1),cex=0.3,bty="n")

#cum base haz:
#ours:
#plot(WB.simData$plot.WB2.H$x,WB.simData$plot.WB2.H$y, type="l", xlim=c(0,c), ylim=c(0,4), ylab=expression(Lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
#orig:
plot(WB.simData$plot.WB2.H$x,WB.simData$plot.WB2.H$y, type="l", xlim=c(0,c), ylim=c(0,2.5), ylab=expression(Lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
for (i in rep){
  lines(t1, y.axis.plot.2.w.H[,i],col="cyan4",lty="dotted")
}
lines(WB.simData$plot.WB2.H, type="l", ylab=expression(lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
lines(t2, mean.2.H,col="red2",lwd=2)
#ours:
# legend(0.1, 3.8, legend=c("truth", "estimates","mean of estimates"),
#        col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.3, bty = "n")
#orig:
legend(0.1, 2.4, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.3, bty = "n")
#3-rd trans:
#base haz:
par(mfrow = c(2,1), mar=c(4.1, 4.1, 1.1, 1.1))
#ours:
# plot(WB.simData$plot.WB3$x,WB.simData$plot.WB3$y, type="l", xlim=c(0,lefttrunc), ylim=c(0,1.1), ylab=expression(lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
#orig:
plot(WB.simData$plot.WB3$x,WB.simData$plot.WB3$y, type="l", xlim=c(0,lefttrunc), ylim=c(0,0.7), ylab=expression(lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
for (i in rep){
  lines(t3, y.axis.plot.3.w.h[,i],col="cyan2",lty="dotted")
}
lines(WB.simData$plot.WB3, type="l", ylab=expression(lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
lines(t3, mean.3.h,col="red2",lwd=2)
#ours:
# legend(0.1, 1, legend=c("truth", "estimates","mean of estimates"),
#        col=c("navy","cyan2","red2"),lty=c(1,3,1),cex=0.32,bty="n")
#orig:
legend(0.1, 0.6, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan2","red2"),lty=c(1,3,1),cex=0.32,bty="n")

#cum base haz:
plot(WB.simData$plot.WB3.H$x,WB.simData$plot.WB3.H$y, type="l", xlim=c(0,lefttrunc), ylim=c(0,1.2), ylab=expression(Lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
for (i in rep){
  lines(t3,y.axis.plot.3.w.H[,i],col="cyan4",lty="dotted")
}
lines(WB.simData$plot.WB3.H, type="l", ylab=expression(lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
lines(t3, mean.3.H,col="red2",lwd=2)
legend(0.1, 1.1, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.32,bty="n")

#run this section if you want to create plots of Bernstein Polynomials:
# BP plots ----------------------------------------------------------------
apply(y.axis.plot.1.b.h,1,mean)->mean.1.BH
apply(y.axis.plot.2.b.h,1,mean)->mean.2.BH
apply(y.axis.plot.3.b.h,1,mean)->mean.3.BH
#1-st trans:
#base haz:
par(mfrow = c(2,1), mar=c(4.1, 4.1, 1.1, 1.1))
plot(WB.simData$plot.WB1$x,WB.simData$plot.WB1$y, type="l", xlim=c(0,c), ylim=c(0,0.13), ylab=expression(lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
rep=1:ncol(y.axis.plot.1.b.h)
for (i in rep){
  lines(t1, y.axis.plot.1.b.h[,i], xlim=c(0,c), ylim=c(0,0.05),col="cyan2",lty="dotted")
}
lines(WB.simData$plot.WB1, type="l", xlim=c(0,c), ylim=c(0,0.07), ylab=expression(lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
lines(t1, mean.1.BH, xlim=c(0,c), ylim=c(0,0.07),col="red2",lwd=2)
legend(0.1, 0.13, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan2","red2"),lty=c(1,3,1),cex=0.3,bty="n")

#cum base haz:
plot(WB.simData$plot.WB1.H$x,WB.simData$plot.WB1.H$y, type="l", xlim=c(0,c), ylim=c(0,2.8), ylab=expression(Lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
for (i in rep){
  h0.interpolate.func <- approxfun(t1, y.axis.plot.1.b.h[,i], rule=2)
  H0 <- sapply(t1, function(x) integrate(h0.interpolate.func, 0, x, stop.on.error = F)$value)
  lines(t1, H0, xlim=c(0,45), ylim=c(0,0.),col="cyan4",lty="dotted")
}
lines(WB.simData$plot.WB1.H$x,WB.simData$plot.WB1.H$y, type="l", xlim=c(0,40), ylim=c(0,1.6), ylab=expression(Lambda[0][1]), xlab=expression(T[1]), col="navy", lwd=2)
h0.interpolate.func <- approxfun(t1, mean.1.BH, rule=2)
H0 <- sapply(t1, function(x) integrate(h0.interpolate.func, 0, x, stop.on.error = F)$value)
lines(t1, H0, xlim=c(0,c), ylim=c(0,1.6),col="red2",lwd=2)
legend(0.1, 2.7, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.3,bty="n")

#2-nd trans:
#base haz:
par(mfrow = c(2,1), mar=c(4.1, 4.1, 1.1, 1.1))
plot(WB.simData$plot.WB2, type="l", xlim=c(0,c), ylim=c(0,0.13), ylab=expression(lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
for (i in rep){
  lines(t2, y.axis.plot.2.b.h[,i], xlim=c(0,30), ylim=c(0,0.25),col="cyan2",lty="dotted")
}
lines(WB.simData$plot.WB2, type="l", xlim=c(0,40), ylim=c(0,0.1), ylab=expression(lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
lines(t2, mean.2.BH, xlim=c(0,40), ylim=c(0,0.08),col="red2",lwd=2)
legend(0.1, 0.13, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan2","red2"),lty=c(1,3,1),cex=0.3,bty="n")

#cum base haz:
plot(WB.simData$plot.WB2.H$x,WB.simData$plot.WB2.H$y, type="l", xlim=c(0,c), ylim=c(0,3.1), ylab=expression(Lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
for (i in rep){
  h0.interpolate.func <- approxfun(t2, y.axis.plot.2.b.h[,i], rule=2)
  H0 <- sapply(t2, function(x) integrate(h0.interpolate.func, 0, x, stop.on.error = F)$value)
  lines(t2, H0, xlim=c(0,40), ylim=c(0,1.1),col="cyan4",lty="dotted")
}
lines(WB.simData$plot.WB2.H$x,WB.simData$plot.WB2.H$y, type="l", xlim=c(0,c), ylim=c(0,1.1), ylab=expression(Lambda[0][2]), xlab=expression(T[2]), col="navy", lwd=2)
h0.interpolate.func <- approxfun(t2, mean.2.BH, rule=2)
H0 <- sapply(t2, function(x) integrate(h0.interpolate.func, 0, x, stop.on.error = F)$value)
lines(t2, H0, xlim=c(0,40), ylim=c(0,0.),col="red2",lwd=2)
legend(0.1, 3.00, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.36,bty="n")

#3-rd trans:
#base haz:
par(mfrow = c(2,1), mar=c(4.1, 4.1, 1.1, 1.1))
plot(WB.simData$plot.WB3$x,WB.simData$plot.WB3$y, type="l", xlim=c(0,7.5), ylim=c(0,0.6), ylab=expression(lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
for (i in rep){
  lines(t3, y.axis.plot.3.b.h[,i], xlim=c(0,15), ylim=c(0,0.5),col="cyan2",lty="dotted")
}
lines(WB.simData$plot.WB3$x,WB.simData$plot.WB3$y, type="l", xlim=c(0,8), ylim=c(0,0.35), ylab=expression(lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
lines(t3, mean.3.BH, xlim=c(0,40), ylim=c(0,0.08),col="red2",lwd=2)
legend(0.1, 0.6, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan2","red2"),lty=c(1,3,1),cex=0.3,bty="n")

#cum base haz:
plot(WB.simData$plot.WB3.H$x,WB.simData$plot.WB3.H$y, type="l", xlim=c(0,7.5), ylim=c(0,1.8), ylab=expression(Lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
for (i in rep){
  h0.interpolate.func <- approxfun(t3, y.axis.plot.3.b.h[,i], rule=2)
  H0 <- sapply(t3, function(x) integrate(h0.interpolate.func, 0, x, stop.on.error = F)$value)
  lines(t3, H0, xlim=c(0,40), ylim=c(0,0.),col="cyan4",lty="dotted")
}
lines(WB.simData$plot.WB3.H$x,WB.simData$plot.WB3.H$y, type="l", xlim=c(0,8), ylim=c(0,0.8), ylab=expression(Lambda[0][3]), xlab=expression(T[2]-T[1]), col="navy", lwd=2)
h0.interpolate.func <- approxfun(t3, mean.3.BH, rule=2)
H0 <- sapply(t3, function(x) integrate(h0.interpolate.func, 0, x, stop.on.error = F)$value)
lines(t3, H0, xlim=c(0,40), ylim=c(0,0.),col="red2",lwd=2)
legend(0.1, 1.8, legend=c("truth", "estimates","mean of estimates"),
       col=c("navy","cyan4","red2"),lty=c(1,3,1),cex=0.3,bty="n")

