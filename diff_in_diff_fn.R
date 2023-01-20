sir.model <- function(t, state, parameters) {
    sir.odes <- function(t,state,parameters) {
        with(as.list(c(state, parameters)),{
            N <- S + I + R
            dS <- -beta.t(t,parameters) * I * S / N
            dI <- beta.t(t,parameters) * I * S / N - gam * I
            dR <- gam * I
            list(c(dS,dI,dR))
        })
    }
    out <- as.data.frame(ode(state,t,sir.odes,parameters))
        r0 <- beta.t(t,parameters)/parameters[["gam"]] 
        out$Rt <- r0 * out$S/(out$S+out$I+out$R)
        out$cumincidence <- out$I + out$R
        
    return(out)
}


beta.t <- function(t,parameters) {
    with(as.list(c(parameters)),{
    if (amp > beta.0) stop("amp should be smaller than beta.0!!!")
    return(beta.0 + amp*sin(t*2*pi/365))
})}


