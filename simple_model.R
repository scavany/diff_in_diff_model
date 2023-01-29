if(!require(pacman)){install.packages('pacman');library(pacman)}
p_load(deSolve,
       adaptivetau)

library(glue)
library(purrr)
library(this.path)

codepath <- this.path::here()

source(glue(codepath, "/diff_in_diff_fn.R"))


## Run models with no intervention -- i.e. we would want our analysis to detect 
 # that there was *no difference* in differences

## Section 1: Same parameters but different initial conditions (same population 
 # size and different proportions infected)

times_1 <- seq(0, 360, 1)
p_1 <- c(
  beta.0 = 1 / 15,
  gam = 1 / 30,
  amp = 0 # i.e. beta is constant over time
)
u0_1a <- c(S = 899, I =  1, R = 100)
u0_1b <- c(S = 890, I = 10, R = 100)

output_1a <- sir.model(times_1, u0_1a, p_1)
output_1b <- sir.model(times_1, u0_1b, p_1)
did_1 <- diffindiff(output_1a, output_1b)
logdid_1 <- logdiffindiff(output_1a, output_1b)

### Difference in differences of point prevalence
plotdid(output_1a$I, output_1b$I, did_1$I, "I")
plotdid(log(output_1a$I), log(output_1b$I), logdid_1$I, "log(I)")

### Difference in differences of R0 
plotdid(output_1a$R0, output_1b$R0, did_1$R0, "R0")

### Difference in differences of Rt
plotdid(output_1a$Rt, output_1b$Rt, did_1$Rt, "Rt")
plotdid(log(output_1a$Rt), log(output_1b$Rt), logdid_1$Rt, "log(Rt)")

### Difference in differences of cumulative incidence 
plotdid(output_1a$cumincidence, output_1b$cumincidence, did_1$cumincidence, "Cumulative incidence")
plotdid(log(output_1a$cumincidence), log(output_1b$cumincidence), logdid_1$cumincidence, "log(Cumulative incidence)")

### Difference in differences of incident cases
plotdid(output_1a$incidentcases, output_1b$incidentcases, did_1$incidentcases, "Incident cases")
plotdid(log(output_1a$incidentcases), log(output_1b$incidentcases), logdid_1$incidentcases, "log(Incident cases)")


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

## Non-time varying parameters. Both pops same parameters. Diff start conditions
t <- seq(from = 0, to = 360, by = .1)
beta <- 1/15
gam <- 1/30
state.1 <- c(S=899,I=1,R=100)
state.2 <- c(S=890,I=10,R=100)
parms.t <- c(beta=beta,gam=gam,beta.0=beta,amp=0/30)

out.det.1 <- sir.model(seq(0,360,0.1),state.1,parms.t)
out.det.2 <- sir.model(seq(0,360,0.1),state.2,parms.t)


pdf("fig2.pdf", pointsize = 16)

par(mfrow=c(1,3))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,xlab="time (days)",ylab="I")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,out.det.1$Rt,type='l',col='red',xlab="Cumulative incidence",ylab="Rt",lwd=3)
lines(out.det.2$R+out.det.2$I,out.det.2$Rt,col='green',lwd=3)

RtvCuminc <- approxfun(out.det.1$R+out.det.1$I,out.det.1$Rt)
plot(out.det.1$R+out.det.1$I,RtvCuminc(out.det.2$R+out.det.2$I)-out.det.2$Rt,xlab="Cumulative incidence",ylab="Difference in Rt estimates")


dev.off()


## Non-time varying parameters. Pops have different parameters. Diff start conditions
t <- seq(from = 0, to = 360, by = .1)
beta <- 1/15
gam <- 1/30
state.1 <- c(S=899,I=1,R=100)
state.2 <- c(S=890,I=10,R=100)
parms.1 <- c(beta=beta,gam=gam,beta.0=1/10,amp=0/30)
parms.2 <- c(beta=beta,gam=gam,beta.0=1/20,amp=0/30)

out.det.1 <- sir.model(seq(0,360,0.1),state.1,parms.1)
out.det.2 <- sir.model(seq(0,360,0.1),state.2,parms.2)


pdf("diffparameters.pdf", pointsize = 16)

par(mfrow=c(1,3))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,xlab="time (days)",ylab="I")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,out.det.1$Rt,type='l',col='red',xlab="Cumulative incidence",ylab="Rt",lwd=3)
lines(out.det.2$R+out.det.2$I,out.det.2$Rt,col='green',lwd=3)

RtvCuminc <- approxfun(out.det.1$R+out.det.1$I,out.det.1$Rt)
plot(out.det.1$R+out.det.1$I,RtvCuminc(out.det.2$R+out.det.2$I)-out.det.2$Rt,xlab="Cumulative incidence",ylab="Difference in Rt estimates")


dev.off()



## Time varying parameters. Both pops same parameters. Diff start conditions
t <- seq(from = 0, to = 360, by = .1)
beta <- 1/15
gam <- 1/30
state.1 <- c(S=899,I=1,R=100)
state.2 <- c(S=890,I=10,R=100)
parms.t <- c(beta=beta,gam=gam,beta.0=beta,amp=1/30)

out.det.1 <- sir.model(seq(0,360,0.1),state.1,parms.t)
out.det.2 <- sir.model(seq(0,360,0.1),state.2,parms.t)


pdf("timevarying.pdf", pointsize = 16)

par(mfrow=c(1,3))
plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,xlab="time (days)",ylab="I")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,out.det.1$Rt,type='l',col='red',xlab="Cumulative incidence",ylab="Rt",lwd=3)
lines(out.det.2$R+out.det.2$I,out.det.2$Rt,col='green',lwd=3)

RtvCuminc <- approxfun(out.det.1$R+out.det.1$I,out.det.1$Rt)
plot(out.det.1$R+out.det.1$I,RtvCuminc(out.det.2$R+out.det.2$I)-out.det.2$Rt,xlab="Cumulative incidence",ylab="Difference in Rt estimates")


dev.off()

lm(outdet.1$RT ~ outdet$cumincidence + t)


##

## Different beta values 
t <- seq(from = 0, to = 360, by = .1)
gamma <- 1/30
state.1 <- c(S=899,I=1,R=100)
state.2 <- c(S=890,I=10,R=100)

parms.t1 <- c(gamma=gamma,beta.0=1/10,amp=1/30)
parms.t2 <- c(gamma=gamma,beta.0=1/20,amp=1/30)

out.det.1 <- sir.model(seq(0,360,0.1),state.1,parms.t1)
out.det.2 <- sir.model(seq(0,360,0.1),state.2,parms.t2)

R0.1 <- beta.t(t, parms.t1)
R0.2 <- beta.t(t, parms.t2)

Rt.1 <- R0.1 * out.det.1$S/(out.det.1$S+out.det.1$I+out.det.1$R)
Rt.2 <- R0.2 * out.det.2$S/(out.det.2$S+out.det.2$I+out.det.2$R)

plot(out.det.1$time,out.det.1$I,col="red",type='l',lwd=3,xlab="time (days)",ylab="I")
lines(out.det.1$time,out.det.2$I,col="green",lwd=3)
plot(out.det.1$R+out.det.1$I,Rt.1,type='l',col='red',xlab="Cumulative incidence",ylab="Rt",lwd=3)
lines(out.det.2$R+out.det.2$I,Rt.2,col='green',lwd=3)
RtvCuminc <- approxfun(out.det.1$R+out.det.1$I,Rt.1)
plot(out.det.1$R+out.det.1$I,RtvCuminc(out.det.2$R+out.det.2$I)-Rt.2,xlab="Cumulative incidence",ylab="Difference in Rt estimates")
