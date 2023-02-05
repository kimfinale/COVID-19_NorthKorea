# particle filter
#
pfilter <- function (params, # parameters
                     y, # initial values of state variables
                     data, # input data set
                     func = NULL, # function to update latent variable 
                     npart = 1000, # number of particles
                     tend = NULL, # simulation stop time
                     dt = 0.2,
                     observed_variable = "daily_symptom",
                     modelled_variable = "CI",
                     backward_sampling = FALSE) {
  
  # Assumptions - using daily growth rate
  nstatevar <- length(y) # number of state variables
  if(is.null(tend)) {
    tend = nrow(data)
  }
  # to store state variables a
  latent_var <- array(0,
                      dim = c(npart, tend, nstatevar),
                      dimnames = list(NULL, NULL, names(y)))
  # latent_var[, 1, ] <- y
  for (nm in names(y)) { # initial value
    latent_var[, 1, nm] <- y[[nm]]
  }
  ## parameters 
  gamma <- params[["gamma"]]
  beta0 <- params[["R0"]] * gamma
  beta_sd <- params[["betavol"]]
  beta <- matrix(rnorm(npart * tend, mean = 0, sd = beta_sd), nrow = tend)
  beta[1,] <- beta0 # this is updated at t=2
  
  wt <- matrix(NA, nrow = npart, ncol = tend) # weight (likelihood)
  wt[, 1] <- 1 / npart  # initial weights
  W <- matrix(NA, nrow = npart, ncol = tend) # normalized weights
  A <- matrix(NA, nrow = npart, ncol = tend) # Resample according to the normalized weight
  
  for (t in 2:tend) {# begin particle loop
    # beta changes according to a Geometric Brownian motion 
    beta[t, ] <- beta[t-1, ] * exp(beta[t, ])
    # run process model
    latent_var[, t, ] <- func(params = params,
                              y = latent_var[, t-1, ],
                              tbegin = t-1,
                              tend = t,
                              dt = dt,
                              beta = beta[t,])
    # calculate weights (likelihood)
    case_expected <- latent_var[, t, modelled_variable]
    case_data <- round(unlist(data[t, observed_variable]))
    expected_val <- pmax(0, case_expected) # make sure that the value is not negative
    # log_lik <- dpois(round(case_data), lambda = expected_val, log = T)
    # wt[, t] <- exp(log_lik)
    wt[, t] <- dpois(round(case_data), lambda = expected_val)
    # normalize particle weights
    W[, t] <- wt[, t] / sum(wt[, t])
    # resample particles by sampling parent particles according to weights
    A[, t] <- sample(1:npart, prob = W[1:npart, t], replace = T)
    # Resample particles for corresponding variables
    latent_var[, t,] <- latent_var[A[, t], t,]
    beta[t,] <- beta[t, A[, t]] #- needed for random walk on beta
  } # end particle loop
  
  # Marginal likelihoods
  lik_values <- rep(NA, tend)
  for (t in 1:tend) {
    lik_values[t] <- log(sum(wt[1:npart, t])) # log-likelihoods
  }# averaged log likelihoods log(L/(npart^tend))
  loglik <- - tend * log(npart) + sum(lik_values)
  
  trace <- replicate(nstatevar + 1, rep(NA, tend), simplify = F)
  names(trace) <- c(names(y), "beta")
  # Latent variables sampled from the filtered distribution (backward sampling)
  if (backward_sampling) {
    loc <- rep(NA, tend)
    # loc_tend <-  sample(1:npart, prob = W[1:npart, tend], replace = T)
    loc[tend] <-  sample.int(npart, size = 1, prob = W[, tend], replace = T)
    
    # particle for the last time step
    trace[["beta"]][tend] <- beta[tend, loc[tend]]
    for (nm in names(y)) {
      trace[[nm]][tend] <- latent_var[loc[tend], tend, nm]
    }
    # lapply(names(y), function(x) trace[[x]][tend] <- latent_var[loc[tend], tend, x])
    # update backward
    for (i in seq(tend, 2, -1)) {
      loc[i-1] <- A[loc[i], i]
      trace[["beta"]][i-1] <- beta[i-1, loc[i-1]]
      
      for (nm in names(y)) {
        trace[[nm]][i-1] <- latent_var[loc[i-1], i-1, nm]
      }
    }
  }
  
  return (list(trace = trace, lik_marginal = lik_values,
               lik_overall_average = loglik,
               latent_var_filtered = latent_var,
               beta_filtered = beta,
               W = W, A = A))
}
