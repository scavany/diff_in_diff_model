if(!require(pacman)){install.packages('pacman');library(pacman)}
p_load(deSolve,
       adaptivetau)

source("./diff_in_diff_fn.R")

state.1 <- c(S=899,I=1,R=100)
state.2 <- c(S=890,I=10,R=100)
beta <- 1/15
gamma <- 1/30
parms <- c(beta=beta,gamma=gamma)
out.det.1 <- sir.model(seq(0,60,0.1),state.1,parms)
out.det.2 <- sir.model(seq(0,60,0.1),state.2,parms)

## Same parameters, different initial conditions
## Difference in difference of prevalence at different time points
## Null intervention at t=0 
par(mfrow=c(1,2))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,log="y")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$time,(out.det.1$I-out.det.2$I)-(state.1[["I"]]-state.2[["I"]]),type='l',lwd=3,xlab="time (days)",ylab="DiD")
abline(h=0)
abline(h=-(state.1[["I"]]-state.2[["I"]]),lty=2)

## Difference in difference of log(prevalence) at different time points
## Null intervention at t=0 
par(mfrow=c(1,2))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,log="y")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
lines(out.det.1$time,out.det.2$I/(state.2[["I"]]/state.1[["I"]]),col="green",lwd=3,lty=2)
plot(out.det.1$time,(log(out.det.1$I)-log(out.det.2$I))-(log(state.1[["I"]])-log(state.2[["I"]])),
     type='l',lwd=3,xlab="time (days)",ylab="DiD")

## Difference in difference of Rt 
par(mfrow=c(1,2))
plot(out.det$time,out.det.1$I,col="red",type='l',lwd=3)
lines(out.det$time,out.det.2$I,col="green",lwd=3)
R0 <- parms[["beta"]]/parms[["gamma"]]
Rt.1 <- R0 * out.det.1$S/(out.det.1$S+out.det.1$I+out.det.1$R)
Rt.2 <- R0 * out.det.2$S/(out.det.2$S+out.det.2$I+out.det.2$R)
plot(out.det.1$time,(Rt.1-Rt.2)-(Rt.1[1]-Rt.2[1]),
     type='l',lwd=3,xlab="time (days)",ylab="DiD")
abline(h=0,lty=2)

par(mfrow=c(1,2))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3)
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,Rt.1,type='l',col='red')
lines(out.det.2$R+out.det.2$I,Rt.2,col='green')

## time-varying beta
parms.t <- c(beta=beta,gamma=gamma,beta.0=beta,amp=1/30)
out.det.1 <- sir.model(seq(0,360,0.1),state.1,parms.t)
out.det.2 <- sir.model(seq(0,360,0.1),state.2,parms.t)
par(mfrow=c(1,2))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3)
lines(out.det.2$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$time,beta.t(out.det.1$time,parms.t[["beta.0"]],parms.t[["amp"]]),type='l',lwd=2)

## Plot the Rt against cumulative incidence
## Check the lines are the same
par(mfrow=c(1,3))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,xlab="time (days)",ylab="I")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,Rt.1,type='l',col='red',xlab="Cumulative incidence",ylab="Rt")
lines(out.det.2$R+out.det.2$I,Rt.2,col='green')
RtvCuminc <- approxfun(out.det.1$R+out.det.1$I,Rt.1)
plot(out.det.1$R+out.det.1$I,RtvCuminc(out.det.2$R+out.det.2$I)-Rt.2,xlab="Cumulative incidence",ylab="Difference in Rt estimates")

## Plot the Rt against cumulative incidence
## Check the lines are the same
parms.t1 <- c(beta=beta,gamma=gamma,beta.0=beta,amp=1/30)
parms.t2 <- c(beta=beta,gamma=gamma,beta.0=beta,amp=1/30)
out.det.1 <- sir.model(seq(0,360,0.1),state.1,parms.t)
out.det.2 <- sir.model(seq(0,360,0.1),state.2,parms.t)
par(mfrow=c(1,3))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,xlab="time (days)",ylab="I")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,Rt.1,type='l',col='red',xlab="Cumulative incidence",ylab="Rt")
lines(out.det.2$R+out.det.2$I,Rt.2,col='green')
RtvCuminc <- approxfun(out.det.1$R+out.det.1$I,Rt.1)
plot(out.det.1$R+out.det.1$I,RtvCuminc(out.det.2$R+out.det.2$I)-Rt.2,xlab="Cumulative incidence",ylab="Difference in Rt estimates")

##

## Different beta values 
parms.t1 <- c(beta=1/10,gamma=gamma,beta.0=1/10,amp=1/30)
parms.t2 <- c(beta=1/20,gamma=gamma,beta.0=1/20,amp=1/30)

out.det.1 <- sir.model(seq(0,360,0.1),state.1,parms.t1)
out.det.2 <- sir.model(seq(0,360,0.1),state.2,parms.t2)


Rt.1 <- R0 * out.det.1$S/(out.det.1$S+out.det.1$I+out.det.1$R)
Rt.2 <- R0 * out.det.2$S/(out.det.2$S+out.det.2$I+out.det.2$R)

plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,xlab="time (days)",ylab="I")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,Rt.1,type='l',col='red',xlab="Cumulative incidence",ylab="Rt")
lines(out.det.2$R+out.det.2$I,Rt.2,col='green')
RtvCuminc <- approxfun(out.det.1$R+out.det.1$I,Rt.1)
plot(out.det.1$R+out.det.1$I,RtvCuminc(out.det.2$R+out.det.2$I)-Rt.2,xlab="Cumulative incidence",ylab="Difference in Rt estimates")
