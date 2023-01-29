sir.model <- function(t, state, parameters) {
  sir.odes <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      N <- S + I + R
      dS <- -beta.t(t, parameters) * I * S / N
      dI <- beta.t(t, parameters) * I * S / N - gam * I
      dR <- gam * I
      list(c(dS, dI, dR))
    })
  }
  out <- as.data.frame(ode(state, t, sir.odes,parameters))
    r0 <- beta.t(t, parameters) / parameters[["gam"]] 
    out$R0 <- r0
    out$Rt <- r0 * out$S/(out$S+out$I+out$R)
    cumincidence <- cuminc(out, state)
    out$cumincidence <- cumincidence
    out$incidentcases <- inc(state, cumincidence)
      
  return(out)
}

cuminc <- function(out, u0) {
  cumincidence <- out$I + out$R - as.list(u0)$R # remove those recovered before start
}

inc <- function(u0, cumincidence) {
  initial = cumincidence[1] - as.list(u0)$I
  incidentcases = c(initial, diff(cumincidence))
  return(incidentcases)
}


## Function for time-varying beta

beta.t <- function(t, parameters) {
  with(as.list(c(parameters)), {
    if (amp > beta.0) stop("amp should be smaller than beta.0!!!")
    return(beta.0 + amp * sin(t * 2 * pi / 365))
  })
}


## Calculate difference in differences between two dataframes

diffindiff <- function(outputa, outputb) {
  diff = outputb - outputa 
  n = nrow(diff)
  sub = diff[rep(2, n), ] # subtract row 2 as row 1 has artefacts due to no new cases
  diffindiff = diff - sub
  return(diffindiff)
}

logdiffindiff <- function(outputa, outputb) {
  diff = log(outputb + 1e-24) - log(outputa + 1e-24) 
  n = nrow(diff)
  sub = diff[rep(2, n), ] # subtract row 2 as row 1 has artefacts due to no new cases
  diffindiff = diff - sub
  return(diffindiff)
}

## Plot values and the difference in differences

plotdid <- function(v1, v2, diff, ylab, ylab2 = "") {
  par(mfrow=c(1, 2))
  plot(v1[-1], col = "red", type = 'l', lwd = 3, xlab = "Time", ylab = ylab)
  lines(v2[-1], col = "blue", type = 'l', lwd = 3)
  plot(diff[-1], type = 'l', lwd = 3, xlab = "Time", ylab = glue("Difference in differences ", ylab2))
  abline(h = 0)
}
