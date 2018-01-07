library(PerformanceAnalytics)
library(dlm)
library(numDeriv)
X<-Design
r<-i
# ols fit - constant parameters
ols.fit = lm(r ~ X)
summary(ols.fit)

# function to build TVP ss model

buildTVP1 <- function(parm, x.mat){
  Design1<-cbind(0,x.mat)
  parm <- exp(parm)
  
  dlm= dlmModReg(X=Design1, addInt=TRUE,
                 dV=parm[2], dW=c(0,parm[3],0,0))
  dlm$GG[1,2]=1
  dlm$GG[2,2]=parm[1]
  
  dlm$m0<-c(2.52045,0,-0.05556,0.56778)
  dlm$C0<-diag(c(3,1,0.05580,0.21317))
  return( dlm )
}

# maximize over log-variances
start.vals = c(-1,0,0)
names(start.vals) = c("lnp", "lns2v", "lns2e")
TVP.mle = dlmMLE(y=r, parm=start.vals,
                 x.mat=X, build=buildTVP1,
                 hessian=T)

dlm
class(TVP.mle)
names(TVP.mle)

TVP.mle
se2 <- sqrt(exp(TVP.mle$par))
names(se2) = c("p", "sv", "se")
sqrt(se2)

#hs<- hessian(function(x) dlmLL(Nile, buildLocalLevel(x)), mleOut$par)
hs<-TVP.mle$hessian
  
all(eigen(hs, only.values = TRUE)$values > 0) # positive definite?

# fitted ss model Filtering with estimated values
TVP.dlm <- buildTVP1(TVP.mle$par, X)
# filtering
TVP.f <- dlmFilter(r, TVP.dlm)
class(TVP.f)

names(TVP.f)

# smoothing
TVP.s <- dlmSmooth(TVP.f)
class(TVP.s)
names(TVP.s)

TVP.s
# extract smoothed states - intercept and slope coefs
alpha.s =TVP.s$s[,1,drop=FALSE]                
beta.s = TVP.s$s[,2,drop=FALSE]
gamma.s=TVP.s$s[,3,drop=FALSE]

colnames(alpha.s) = "alpha"
colnames(beta.s) = "beta"
colnames(gamma.s) = "gamma"

# extract std errors - dlmSvd2var gives list of MSE matrices
  mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)

se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
se.xts = ts(se.mat,  frequency=12, start=c(2002,1))
colnames(se.xts) = c("alpha", "beta","","")
a.u = alpha.s + 1.96*se.xts[,"alpha"]
a.l = alpha.s - 1.96*se.xts[, "alpha"]
b.u = beta.s + 1.96*se.xts[,"beta"]
b.l = beta.s - 1.96*se.xts[, "beta"]
c.u = gamma.s + 1.96*se.xts[,"gamma"]
c.l = gamma.s - 1.96*se.xts[, "gamma"]

e=alpha.s- r
e.u = e + 1.96*se.xts[,"alpha"]
e.l = e - 1.96*se.xts[, "alpha"]
chart.TimeSeries(cbind(alpha.s, a.l, a.u), main="Smoothed estimates of alpha",
                colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(alpha),xlab="")
fitted<-c(0.07390197, -0.04775924)*X
react<-rowSums(fitted)
react1<--0.04775924*X[,2]
react2<-0.07390197*X[,1]
chart.TimeSeries(cbind(e,react2), main="Smoothed estimates of alpha and observed values",
                  colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(alpha),xlab="")

chart.TimeSeries(cbind(beta.s, b.l, b.u), main="Smoothed estimates of beta",
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(beta),xlab="")
chart.TimeSeries(cbind(gamma.s, c.l, c.u), main="Smoothed estimates of gamma",
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(beta),xlab="")


library(forecast)
tsdisplay(diff(r),main="")
Design1<-cbind(0,Design)
tvp.dlm = dlmModReg(X=Design1, addInt=TRUE,
                    dV=s2v, dW=c(0,s2a,0,0))
tvp.dlm$G[1,2]=1
tvp.dlm$G[2,2]=p


buildTVP1 <- function(parm, x.mat){
  Design1<-cbind(0,x.mat)
  parm <- exp(parm)
  
  dlm= dlmModReg(X=Design1, addInt=TRUE,
                      dV=parm[2], dW=c(0,parm[3],0,0))
  dlm$GG[1,2]=1
  dlm$GG[2,2]=parm[1]
  
  return( dlm )
}

# with one regressor
#####################################################
#####
ols.fit = lm(r ~ X[,2])
summary(ols.fit)


buildTVP2 <- function(parm, x.mat){
  Design1<-x.mat
  parm <- exp(parm)
  
  dlm= dlmModReg(X=Design, addInt=TRUE,
                 dV=parm[2], dW=c(parm[3],0,0))
  dlm$GG[1,1]=parm[1]
  
  
  dlm$m0<-c(2.52045,0.07,0.56778)
  dlm$C0<-diag(c(3,1,0.21317))
  return( dlm )
}
start.vals = c(-1,0,0)
names(start.vals) = c("lnp", "lns2v", "lns2e")
TVP.mle = dlmMLE(y=r, parm=start.vals,
                 x.mat=X, build=buildTVP2,
                 hessian=T)


TVP.mle
se2 <- exp(TVP.mle$par)
names(se2) = c("p", "sv", "se")

se2
#hs<- hessian(function(x) dlmLL(Nile, buildLocalLevel(x)), mleOut$par)
hs<-TVP.mle$hessian

all(eigen(hs, only.values = TRUE)$values > 0) # positive definite?
aVar <- solve(hs) # asymptotic variance/covariance matrix
sqrt(diag(aVar))

# fitted ss model Filtering with estimated values
TVP.dlm <- buildTVP2(TVP.mle$par, X)
# filtering
TVP.f <- dlmFilter(r, TVP.dlm)

# smoothing
TVP.s <- dlmSmooth(TVP.f)


TVP.s$s
# extract smoothed states - intercept and slope coefs
alpha.s =TVP.s$s[,1,drop=FALSE]                
beta.s = TVP.s$s[,2,drop=FALSE]
gamma.s=TVP.s$s[,3,drop=FALSE]

colnames(alpha.s) = "alpha"
colnames(beta.s) = "beta"
colnames(gamma.s) = "gamma"

> # extract std errors - dlmSvd2var gives list of MSE matrices
  mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)

se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
se.xts = ts(se.mat,  frequency=12, start=c(2002,1))
colnames(se.xts) = c("alpha", "beta","","")
a.u = alpha.s + 1.96*se.xts[,"alpha"]
a.l = alpha.s - 1.96*se.xts[, "alpha"]
b.u = beta.s + 1.96*se.xts[,"beta"]
b.l = beta.s - 1.96*se.xts[, "beta"]
c.u = gamma.s + 1.96*se.xts[,"gamma"]
c.l = gamma.s - 1.96*se.xts[, "gamma"]

e=alpha.s- r
e.u = e + 1.96*se.xts[,"alpha"]
e.l = e - 1.96*se.xts[, "alpha"]
chart.TimeSeries(cbind(alpha.s, a.l, a.u), main="Smoothed estimates of alpha",
                 , colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(alpha),xlab="")

chart.TimeSeries(cbind(e), main="Smoothed estimates of alpha and observed values",
                 , colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(alpha),xlab="")

chart.TimeSeries(cbind(beta.s, b.l, b.u), main="Smoothed estimates of beta",
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(beta),xlab="")
chart.TimeSeries(cbind(gamma.s-e4, c.l-e4, c.u-e4), main="Smoothed estimates of gamma",
                 colorset=c(1,2,2), lty=c(1,2,2),ylab=expression(beta),xlab="")

