stdDev <- 1;
 x <- seq(-5,5,by=0.01)
 y <- dnorm(x,sd=stdDev)
 right <- qnorm(0.05,sd=stdDev)
 plot(x,y,type="l",xaxt="n",ylab="Funds",
       xlab=expression(paste('Assumed Distribution of investment oppourtunities')),
       axes=FALSE,ylim=c(0,max(y)*1.05),xlim=c(min(x),max(x)),
       frame.plot=FALSE)
 axis(1,at=c(-5,right,0,5),
       pos = c(0,0),
       labels=c(expression(' '),expression(bar(i)[policy]),expression(mu[0]),expression(' ')))
 axis(2)
 xReject <- seq(-5,right,by=0.01)
 yReject <- dnorm(xReject,sd=stdDev)
 
polygon(c(xReject,xReject[length(xReject)],xReject[1]),
          c(yReject,0, 0), col='red')

polygon(c(xReject,xReject[length(xReject)],xReject[1]),
        c(yReject,0, 0), col='red')
legend("topright", c(" Funds 'invested' in the Central Bank"),fill="red")





########################################################################
#Impuls responses from one dimensional VECM.

IRF<-function(Parm,n,y0){
 beta<-Parm[1]
 EC<-Parm[2]
 alfa1<-Parm[3]
 alfa2<-Parm[4]
 const<-Parm[5]
  x<-rep(0,n)
 y<-rep(y0,n)
 y[2:n]=y0+0.01*abs(y0)
 dx<-rep(0,n)
 dy<-rep(0,n)
 dy[2]=y[2]-y[1]
 
 x[1]=+beta*y[1]-const/EC
 for(i in 2:n){
    dx[i]=EC*(x[i-1]-beta*y[i-1])+alfa1*dy[i]+alfa2*dx[i-1]+const
    x[i]=x[i-1]+dx[i]
  }
  return(100*(x/x[1]-1))
}
IRF1<-function(Parm,n,y0){
  beta<-Parm[1]
  EC<-Parm[2]
  alfa1<-Parm[3]
  alfa2<-Parm[4]
  const<-Parm[5]
  x<-rep(0,n)
  y<-rep(y0,n)
  y[2:n]=y0+0.01
  dx<-rep(0,n)
  dy<-rep(0,n)
  dy[2]=y[2]-y[1]
  
  x[1]=+beta*y[1]-const/EC
  for(i in 2:n){
    dx[i]=EC*(x[i-1]-beta*y[i-1])+alfa1*dy[i]+alfa2*dx[i-1]+const
    x[i]=x[i-1]+dx[i]
  }
  return(100*(x-x[1])[2:n])
}

#Real:
Parm<-rep(0.001,5)
Parm1<-c(0.510,-0.07,0.492,0.144,0.116)
Parm2<-c(0.419,-0.07,0.498,0.145,0.122)
Parm3<-c(0.551,-0.10,0.812,0.150,-0.23)
Parm4<-c(0.769,-0.10,0.499,0.165,0.139)



ts.plot(ts(cbind(IRF1(Parm1,36,1),
                 IRF1(Parm2,36,1),
                 IRF1(Parm3,36,1),
                 IRF1(Parm4,36,1))),gpars= list(col=1:4),main="Impulse response functions for change of the European Real Rates")
legend("topright",c("PMG","MG","AMG","DFE") , col=1:4, lty=1, cex=.65)



#Nominal
Parm1<-c(2.055,-0.01,0.532,0.301,0.050)
Parm2<-c(0.857,-0.03,0.528,0.296,0.029)
Parm3<-c(3.121,-0.11,1.002,0.300,-0.03)
Parm4<-c(1.654,-0.05,0.543,0.228,0.048)

#bez AMG
ts.plot(ts(cbind(IRF1(Parm1,36,1),
                 IRF1(Parm2,36,1),
                 IRF1(Parm4,36,1))),gpars= list(col=c(1,3,4)), main="Impulse response functions- nominal Rates" )
legend("topright",c("PMG","MG","DFE") , col=c(1,3,4), lty=1, cex=.65)

