# Update the SEPAIR state variables over the time step dt
SEPAIR_step <- function (params = NULL,
                       y = NULL,
                       tbegin = 0,
                       tend = 1,
                       dt = 0.2,
                       beta = NULL) {
  
  S <- y[, "S"]
  E <- y[, "E"]
  P <- y[, "P"]
  A <- y[, "A"]
  I <- y[, "I"]
  R <- y[, "R"]
  
  daily_infected <- rep(0, nrow(y)) # set to zero to store the change over the dt
  daily_symptom <- rep(0, nrow(y))
  
  
  ## set beta based on the predefined Rt
  ## first set the duration of infectiousness correct
  ## account for the relative infectiousness of P and A states
  # durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]])
  # durI <- (1 / params[["gamma"]])
  fa <- params[["fa"]]# fraction of asymptomatic state
  bp <- params[["bp"]] # relative infectiousness of pre-symptomatic state
  ba <- params[["ba"]] # relative infectiousness of asymptomatic state
  # R0_dur <- ((1 - fa) + fa * ba) * durI + bp * durP
  
  ## extract the parameters
  epsilon <- params[["epsilon"]]
  delta <- params[["delta"]]
  gamma <- params[["gamma"]]
  
  N <- S + E + P + A + I + R
  
  for (i in seq((tbegin + dt), tend, dt)) {
    # beta is already assumed to be adjusted by N such that it can
    # be translated to Rt by multiplying the duration of infectiousness
    
    rate_from_p <- 1 / (1/delta - 1/epsilon) # rate of transitioning from state P
    # change of state variables over dt
    StoE <- beta * (bp * P + ba * A + I) * dt
    EtoP <- epsilon * E * dt
    PtoA <- rate_from_p * fa * P * dt
    PtoI <- rate_from_p * (1 - fa) * P * dt
    AtoR <- gamma * A * dt
    ItoR <- gamma * I * dt
    
    # Update state variables
    S <- S - StoE
    E <- E + StoE - EtoP
    P <- P + EtoP - PtoA - PtoI
    A <- A + PtoA - AtoR
    I <- I + PtoI - ItoR
    R <- R + AtoR + ItoR
    
    daily_infected <- daily_infected + StoE
    daily_symptom <- daily_symptom + PtoI
  }
  
  y[, "S"] <- S
  y[, "E"] <- E
  y[, "P"] <- P
  y[, "A"] <- A
  y[, "I"] <- I
  y[, "R"] <- R
  y[, "CE"] <- daily_infected
  y[, "CI"] <- daily_symptom
  
  return(y)
}
