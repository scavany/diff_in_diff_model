sir.model <- function(t, state, parameters)
{
    sir.odes <- function(t,state,parameters) {
        with(as.list(c(state, parameters)),{
            N <- S + I + R
            dS <- -beta.t(t,beta.0,amp) * I * S / N
            dI <- beta.t(t,beta.0,amp) * I * S / N - gamma * I
            dR <- gamma * I
            list(c(dS,dI,dR))
        })
    }
    out <- as.data.frame(ode(state,t,sir.odes,parameters))
    return(out)
}

beta.t <- function(t,beta.0,amp) {
    if (amp > beta.0) stop("amp should be smaller than beta.0!!!")
    return(beta.0 + amp*sin(t*2*pi/365))
}
