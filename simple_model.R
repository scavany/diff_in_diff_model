if(!require(pacman)){install.packages('pacman');library(pacman)}
p_load(deSolve,
       adaptivetau)

sir.model <- function(t, state, parameters)
{
    sir.odes <- function(t,state,parameters) {
        with(as.list(c(state, parameters)),{
            N <- S + I + R
            dS <- -beta * I * S / N
            dI <- beta * I * S / N - gamma * I
            dR <- gamma * I
            list(c(dS,dI,dR))
        })
    }
    out <- as.data.frame(ode(state,t,sir.odes,parameters))
    return(out)
}

sir.stoch <- function(t, state, parameters){
    ## Transitions matrix
    sir.transitions <- list(
        c(I=1,S=-1),
        c(R=1,I=-1)
    )
    ## Rates function
    sir.rates <- function(state, parameters, t)
    {
        with(as.list(c(state, parameters)),{
            N <- S + I + R
            return(c(beta * I * S / N,
                     gamma * I))
        })
    }
    
    out.temp <- as.data.frame(adaptivetau::ssa.adaptivetau(state,  sir.transitions,
                                                           sir.rates,
                                                           parameters, tf=diff(range(t))))
    out.temp$time <- out.temp$time + min(t)
    out <- cbind(time = t, apply(out.temp[, -1], 2, function(col) {
        approx(x = out.temp[, 1], y = col, xout = t, method = "constant")$y
    }))
    return(as.data.frame(out))
}


state.init <- c(S=9999,I=1,R=0)
beta <- 1
gamma <- 1/5
parms <- c(beta=beta,gamma=gamma)
out.det <- sir.model(seq(0,60,0.1),state.init,parms)
out.stoch <- sir.stoch(seq(0,60,0.1),state.init,parms)

plot(out.det$time,out.det$S,col="red",type='l',lwd=3)
lines(out.det$time,out.det$I,col="green",lwd=3)
lines(out.det$time,out.det$R,col="blue",lwd=3)
lines(out.stoch$time,out.stoch$S,col="red",lwd=3,lty=2)
lines(out.stoch$time,out.stoch$I,col="green",lwd=3,lty=2)
lines(out.stoch$time,out.stoch$R,col="blue",lwd=3,lty=2)
