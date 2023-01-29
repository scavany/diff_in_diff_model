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

## Scenario 1: Same parameters but different initial conditions (same population 
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
plotdid(output_1a$cumincidence, output_1b$cumincidence, did_1$cumincidence, 
  "Cumulative incidence")
plotdid(log(output_1a$cumincidence), log(output_1b$cumincidence), logdid_1$cumincidence, 
  "log(Cumulative incidence)")

### Difference in differences of incident cases
plotdid(output_1a$incidentcases, output_1b$incidentcases, did_1$incidentcases, 
  "Incident cases")
plotdid(log(output_1a$incidentcases), log(output_1b$incidentcases), logdid_1$incidentcases, 
  "log(Incident cases)")


## Scenario 2: Same parameters but different initial conditions (different population 
# size and same proportions initially infected)

times_2 <- seq(0, 360, 1)
p_2 <- c(
  beta.0 = 1 / 15,
  gam = 1 / 30,
  amp = 0 # i.e. beta is constant over time
)
u0_2a <- c(S =  898, I = 2, R = 100)
u0_2b <- c(S = 1347, I = 3, R = 150)

output_2a <- sir.model(times_2, u0_2a, p_2)
output_2b <- sir.model(times_2, u0_2b, p_2)
did_2 <- diffindiff(output_2a, output_2b)
logdid_2 <- logdiffindiff(output_2a, output_2b)

### Difference in differences of point prevalence
plotdid(output_2a$I, output_2b$I, did_2$I, "I")
plotdid(log(output_2a$I), log(output_2b$I), logdid_2$I, "log(I)")

### Difference in differences of R0 
plotdid(output_2a$R0, output_2b$R0, did_2$R0, "R0")

### Difference in differences of Rt
plotdid(output_2a$Rt, output_2b$Rt, did_2$Rt, "Rt")
plotdid(log(output_2a$Rt), log(output_2b$Rt), logdid_2$Rt, "log(Rt)")

### Difference in differences of cumulative incidence 
plotdid(output_2a$cumincidence, output_2b$cumincidence, did_2$cumincidence, 
  "Cumulative incidence")
plotdid(log(output_2a$cumincidence), log(output_2b$cumincidence), logdid_2$cumincidence, 
  "log(Cumulative incidence)")

### Difference in differences of incident cases
plotdid(output_2a$incidentcases, output_2b$incidentcases, did_2$incidentcases, 
  "Incident cases")
plotdid(log(output_2a$incidentcases), log(output_2b$incidentcases), logdid_2$incidentcases, 
  "log(Incident cases)")


## Scenario 3: Same parameters, time-varying beta, but different initial conditions 
# (same population size and different proportions infected)

times_3 <- seq(0, 360, 1)
p_3 <- c(
  beta.0 = 1 / 15,
  gam = 1 / 30,
  amp = .3 / 15
)
u0_3a <- c(S = 899, I =  1, R = 100)
u0_3b <- c(S = 890, I = 10, R = 100)

output_3a <- sir.model(times_3, u0_3a, p_3)
output_3b <- sir.model(times_3, u0_3b, p_3)
did_3 <- diffindiff(output_3a, output_3b)
logdid_3 <- logdiffindiff(output_3a, output_3b)

### Difference in differences of point prevalence
plotdid(output_3a$I, output_3b$I, did_3$I, "I")
plotdid(log(output_3a$I), log(output_3b$I), logdid_3$I, "log(I)")

### Difference in differences of R0 
plotdid(output_3a$R0, output_3b$R0, did_3$R0, "R0")

### Difference in differences of Rt
plotdid(output_3a$Rt, output_3b$Rt, did_3$Rt, "Rt")
plotdid(log(output_3a$Rt), log(output_3b$Rt), logdid_3$Rt, "log(Rt)")

### Difference in differences of cumulative incidence 
plotdid(output_3a$cumincidence, output_3b$cumincidence, did_3$cumincidence, 
        "Cumulative incidence")
plotdid(log(output_3a$cumincidence), log(output_3b$cumincidence), logdid_3$cumincidence, 
        "log(Cumulative incidence)")

### Difference in differences of incident cases
plotdid(output_3a$incidentcases, output_3b$incidentcases, did_3$incidentcases, 
        "Incident cases")
plotdid(log(output_3a$incidentcases), log(output_3b$incidentcases), logdid_3$incidentcases, 
        "log(Incident cases)")


## Scenario 4: Same parameters, time-varying betas, but different initial conditions 
 # (different population size and same proportions initially infected)

times_4 <- seq(0, 360, 1)
p_4 <- c(
  beta.0 = 1 / 15,
  gam = 1 / 30,
  amp = .3 / 15
)
u0_4a <- c(S =  898, I = 2, R = 100)
u0_4b <- c(S = 1347, I = 3, R = 150)

output_4a <- sir.model(times_4, u0_4a, p_4)
output_4b <- sir.model(times_4, u0_4b, p_4)
did_4 <- diffindiff(output_4a, output_4b)
logdid_4 <- logdiffindiff(output_4a, output_4b)

### Difference in differences of point prevalence
plotdid(output_4a$I, output_4b$I, did_4$I, "I")
plotdid(log(output_4a$I), log(output_4b$I), logdid_4$I, "log(I)")

### Difference in differences of R0 
plotdid(output_4a$R0, output_4b$R0, did_4$R0, "R0")

### Difference in differences of Rt
plotdid(output_4a$Rt, output_4b$Rt, did_4$Rt, "Rt")
plotdid(log(output_4a$Rt), log(output_4b$Rt), logdid_4$Rt, "log(Rt)")

### Difference in differences of cumulative incidence 
plotdid(output_4a$cumincidence, output_4b$cumincidence, did_4$cumincidence, 
        "Cumulative incidence")
plotdid(log(output_4a$cumincidence), log(output_4b$cumincidence), logdid_4$cumincidence, 
        "log(Cumulative incidence)")

### Difference in differences of incident cases
plotdid(output_4a$incidentcases, output_4b$incidentcases, did_4$incidentcases, 
        "Incident cases")
plotdid(log(output_4a$incidentcases), log(output_4b$incidentcases), logdid_4$incidentcases, 
        "log(Incident cases)")


## Scenario 5: Degree of time-varying beta differs but other parameters the same,
 # same initial conditions 

times_5 <- seq(0, 360, 1)
p_5a <- c(
  beta.0 = 1 / 15,
  gam = 1 / 30,
  amp = .3 / 15
)
p_5b <- c(
  beta.0 = 1 / 15,
  gam = 1 / 30,
  amp = .6 / 15
)
u0_5 <- c(S = 899, I =  1, R = 100)

output_5a <- sir.model(times_5, u0_5, p_5a)
output_5b <- sir.model(times_5, u0_5, p_5b)
did_5 <- diffindiff(output_5a, output_5b)
logdid_5 <- logdiffindiff(output_5a, output_5b)

### Difference in differences of point prevalence
plotdid(output_5a$I, output_5b$I, did_5$I, "I")
plotdid(log(output_5a$I), log(output_5b$I), logdid_5$I, "log(I)")

### Difference in differences of R0 
plotdid(output_5a$R0, output_5b$R0, did_5$R0, "R0")

### Difference in differences of Rt
plotdid(output_5a$Rt, output_5b$Rt, did_5$Rt, "Rt")
plotdid(log(output_5a$Rt), log(output_5b$Rt), logdid_5$Rt, "log(Rt)")

### Difference in differences of cumulative incidence 
plotdid(output_5a$cumincidence, output_5b$cumincidence, did_5$cumincidence, 
        "Cumulative incidence")
plotdid(log(output_5a$cumincidence), log(output_5b$cumincidence), logdid_5$cumincidence, 
        "log(Cumulative incidence)")

### Difference in differences of incident cases
plotdid(output_5a$incidentcases, output_5b$incidentcases, did_5$incidentcases, 
        "Incident cases")
plotdid(log(output_5a$incidentcases), log(output_5b$incidentcases), logdid_5$incidentcases, 
        "log(Incident cases)")


## Scenario 6: Beta parameters differ (not time-varying) but other parameters the same,
# same initial conditions 

times_6 <- seq(0, 360, 1)
p_6a <- c(
  beta.0 = 1 / 15,
  gam = 1 / 30,
  amp = 0
)
p_6b <- c(
  beta.0 = 1.3 / 15,
  gam = 1 / 30,
  amp = 0
)
u0_6 <- c(S =  898, I = 2, R = 100)

output_6a <- sir.model(times_6, u0_6, p_6a)
output_6b <- sir.model(times_6, u0_6, p_6b)
did_6 <- diffindiff(output_6a, output_6b)
logdid_6 <- logdiffindiff(output_6a, output_6b)

### Difference in differences of point prevalence
plotdid(output_6a$I, output_6b$I, did_6$I, "I")
plotdid(log(output_6a$I), log(output_6b$I), logdid_6$I, "log(I)")

### Difference in differences of R0 
plotdid(output_6a$R0, output_6b$R0, did_6$R0, "R0")

### Difference in differences of Rt
plotdid(output_6a$Rt, output_6b$Rt, did_6$Rt, "Rt")
plotdid(log(output_6a$Rt), log(output_6b$Rt), logdid_6$Rt, "log(Rt)")

### Difference in differences of cumulative incidence 
plotdid(output_6a$cumincidence, output_6b$cumincidence, did_6$cumincidence, 
        "Cumulative incidence")
plotdid(log(output_6a$cumincidence), log(output_6b$cumincidence), logdid_6$cumincidence, 
        "log(Cumulative incidence)")

### Difference in differences of incident cases
plotdid(output_6a$incidentcases, output_6b$incidentcases, did_6$incidentcases, 
        "Incident cases")
plotdid(log(output_6a$incidentcases), log(output_6b$incidentcases), logdid_6$incidentcases, 
        "log(Incident cases)")






