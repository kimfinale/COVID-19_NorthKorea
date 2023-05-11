# extract_trace
extract_trace <- function (params = NULL, 
                           y = NULL, 
                           data = NULL, 
                           npart = 1000, 
                           tend = 200, 
                           dt = 0.2,
                           func = NULL,
                           observed_variable = "daily_infected", 
                           modelled_variable = "CE",
                           backward_sampling = TRUE) {
  nstatevar <- length(y)
  if (is.data.frame(y)) {
    nstatevar <- ncol(y)
  }
  
  res <- pfilter(params = params,
          y = y, # initial values of state variables
          func = func,
          npart = npart,
          data = data,
          tend = tend,
          dt = dt,
          observed_variable = observed_variable, 
          modelled_variable = modelled_variable,
          backward_sampling = backward_sampling) # input data set)
  
  
  output <- data.frame(matrix(NA, nrow = tend, ncol = nstatevar + 
                                1))
  names(output) <- c(names(y), "Rt")
  for (nm in names(y)) {
    output[, nm] <- res$trace[[nm]]
  }
  # output[, "beta"] <- res$trace[["beta"]]
  
  durP <- (1/params[["delta"]] - 1/params[["epsilon"]])
  durI <- (1/params[["gamma"]])
  fa <- params[["fa"]]
  ba <- params[["ba"]]
  bp <- params[["bp"]]
  R0_dur <- ((1 - fa) + fa * ba) * durI + bp * durP
  output[, "Rt"] <- res$trace[["beta"]] * R0_dur
  return(output)
}
