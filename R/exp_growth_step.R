# Update the SEPAIR state variables over the time step dt
exp_growth_step <- function (params = NULL,
                      y0 = NULL,
                      tbegin = 0,
                      tend = 100,
                      dt = 0.1,
                      saveat = 1) {

  y = data.frame(matrix(NA, nrow=ceiling((tend-tbegin)/saveat)+1,
                        ncol=length(y0)+1))
  names(y) = c("t", names(y0))
  y$t <- seq(0, length.out=nrow(y), by=saveat)
  y[1, names(y0)] = y0

  r <- params[["r"]]

  I = y0[[1]]

  for (i in 1:nrow(y)) {
    for (j in seq(dt, saveat, dt)) {
      I <- I + I * r * dt
    }
    y[i, "I"] <- I
  }

  return(y)
}
