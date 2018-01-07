library(PerformanceAnalytics)
library(foreign)
library(dlm)
data<-read.dta("wig20ver3.dta") 
Interwencja <- read.csv("Interwencja.csv", sep=";")
data<-cbind(data,100* Interwencja[2],100*Interwencja[3],100*Interwencja[4], data$interestrate-data$inflation)
names(data)[36]<-"realinterest"
tsdata1<-ts(data[2:36],  frequency=12, start=c(2000,1))

# real interest rate
#r=diff(tsdata1[,c("realinterest")])
r=tsdata1[,c("realinterest")]
r<-r[25:189]
r<-ts(r,  frequency=12, start=c(2002,1))
i<-tsdata1[,c("interestrate")]
i<-i[25:189]
i<-ts(i,  frequency=12, start=c(2002,1))

pi<-tsdata1[,c("inflation")]
pi<-pi[25:189]
pi<-ts(pi,  frequency=12, start=c(2002,1))

# design Matrix
X=tsdata1[,c("Interwencja")]
X=cbind(X-mean(X),diff(X))
Design<-ts(X[25:189,],  frequency=12, start=c(2002,1))
colnames(Design)<-c("Interwencja","dInterwencja")

 
s2v = 1
 s2a = 0.01
 s2b = 0.01
s2c = 0.01
 tvp.dlm = dlmModReg(X=Design, addInt=TRUE,
                      dV=s2v, dW=c(s2a, s2b,s2c))

tvp.dlm[c("FF","V","GG","W","m0","C0")]
tvp.dlm[c("JFF","JV","JGG","JW")]
head(tvp.dlm$X)

# ols fit - constant parameters
 ols.fit = lm(r ~ Design)
 summary(ols.fit)

 # function to build TVP ss model
 buildTVP <- function(parm, x.mat){
    parm <- exp(parm)
    return( dlmModReg(X=x.mat, dV=parm[1],
                      dW=c(parm[2], 0,0)) )
  }
 # maximize over log-variances
 start.vals = c(0,0)

 TVP.mle = dlmMLE(y=r, parm=start.vals,
                   x.mat=Design, build=buildTVP,
                   hessian=T)

class(TVP.mle)
names(TVP.mle)

TVP.mle
se2 <- sqrt(exp(TVP.mle$par))
 names(se2) = c("sv", "sa", "sb","sc")
 se2

# fitted ss model Filtering with estimated values
 TVP.dlm <- buildTVP(TVP.mle$par, Design)
 # filtering
 TVP.f <- dlmFilter(r, TVP.dlm)
 class(TVP.f)

names(TVP.f)

# smoothing
TVP.s <- dlmSmooth(TVP.f)
 class(TVP.s)
names(TVP.s)


# extract smoothed states - intercept and slope coefs
 alpha.s =TVP.s$s[,1,drop=FALSE]                
 beta.s = TVP.s$s[,2,drop=FALSE]
  gamma.s=TVP.s$s[,3,drop=FALSE]
TVP.s$s          
colnames(alpha.s) = "alpha"
 colnames(beta.s) = "beta"
colnames(gamma.s) = "gamma"

 # extract std errors - dlmSvd2var gives list of MSE matrices
   mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)
 se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
 se.xts = ts(se.mat,  frequency=12, start=c(2002,1))
 colnames(se.xts) = c("alpha", "beta","gamma")
 a.u = alpha.s + 1.96*se.xts[,"alpha"]
 a.l = alpha.s - 1.96*se.xts[, "alpha"]
 b.u = beta.s + 1.96*se.xts[,"beta"]
 b.l = beta.s - 1.96*se.xts[, "beta"]
c.u = gamma.s + 1.96*se.xts[,"gamma"]
c.l = gamma.s - 1.96*se.xts[, "gamma"]

e=alpha.s-r
chart.TimeSeries(cbind(alpha.s, a.l, a.u), main="Smoothed estimates of alpha",
                  colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(alpha),xlab="")

chart.TimeSeries(cbind(alpha.s,r), main="Smoothed estimates of alpha and observed values",
                  col=c("red","green"), lty=c(1,2),ylab=expression(alpha),xlab="")
chart.TimeSeries(cbind(e), main="Smoothed estimates of changed interest rate",
                  col=c("red","green"), lty=c(1,2),ylab=expression(alpha),xlab="")
legend("topright", c("Real interest rate(Observed)", "Natural real interest rate (filtered)"),col=c("red","green"),  lty=1, cex=.65)
chart.TimeSeries(cbind(beta.s, b.l, b.u), main="Smoothed estimates of beta",
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(beta),xlab="")
chart.TimeSeries(cbind(gamma.s, c.l, c.u), main="Smoothed estimates of gamma",
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(beta),xlab="")
