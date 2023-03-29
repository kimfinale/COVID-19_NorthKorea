# Update the SEPAIR state variables over the time step dt
SIR_step <- function (params = NULL,
                      y0 = NULL,
                      tbegin = 0,
                      tend = 100,
                      dt = 0.1,
                      saveat = 1) {

  y = data.frame(matrix(NA, nrow=ceiling((tend-tbegin)/saveat)+1,
                        ncol=(length(y0)+2)))
  names(y) = c("t", names(y0), "CI")
  y$t <- seq(0, length.out=nrow(y), by=saveat)
  y[1, names(y0)] = y0
  y[1, "CI"] = 0

  cumul_infected <- 0

  gamma <- params[["gamma"]]
  beta <- params[["beta"]]

  S = y0[[1]]
  I = y0[[2]]
  R = y0[[3]]

  for (i in 1:nrow(y)) {
    for (j in seq(dt, saveat, dt)) {
      N <- S + I + R
      StoI <- beta * S * I/N * dt
      ItoR <- gamma * I * dt

      # Update state variables
      S <- S - StoI
      I <- I + StoI - ItoR
      R <- R + ItoR

      cumul_infected <- cumul_infected + StoI
    }

    y[i, "S"] <- S
    y[i, "I"] <- I
    y[i, "R"] <- R
    y[i, "CI"] <- cumul_infected
  }

  return(y)
}
