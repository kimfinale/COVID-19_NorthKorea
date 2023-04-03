#' This function implements forward filtering-backward sampling algorithm
#'
#' The \code{particle_filter()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param params Parameters
#' @param y A vector of state variables (e.g., c(S,E,P,I,R))
#' @param data data to calculate the likelihoods against
#' @param npart number of particles
#' @export
#' @import data.table
#' @examples
#' sample <- particle_filter()
pfilter <- function (params = theta,
                     y = y0,
                     data = NULL,
                     npart = 1000,
                     tend = 200,
                     backward_sampling = TRUE) {

  # Assumptions - using daily growth rate
  nstatevar <- length(y)
  latent_var <- array(0,
                      dim = c(npart, tend, nstatevar),
                      dimnames = list(NULL, NULL, names(y)))
  # latent_var[, 1, ] <- y
  for (nm in names(y)) {
    latent_var[, 1, nm] <- y[[nm]]
  }

  beta_sd <- params[["betavol"]]
  beta <- matrix(rnorm(npart * tend, mean = 0, sd = beta_sd), nrow = tend)
  beta[1,] <- beta0 # this is updated at t=2

  wt <- matrix(NA, nrow = npart, ncol = tend) # weight (likelihood)
  wt[, 1] <- 1 / npart  # initial weights
  W <- matrix(NA, nrow = npart, ncol = tend) # normalized weights
  A <- matrix(NA, nrow = npart, ncol = tend) # Resample according to the normalized weight

  for (t in 2:tend) {# begin particle loop
    # cat("t =", t)
    # geometric
    beta[t, ] <- beta[t-1, ] * exp(beta[t, ])
    # cat(", sd =", beta_sd, "\n")
    # run process model
    latent_var[, t, ] <- process_model(y=latent_var[, t-1, ],
                                       beta = beta[t,])
    # calculate weights (likelihood)
    wt[, t] <- assign_weights(var=latent_var, t=t, data=data)
    # normalize particle weights
    # cat("sum(is.na(wt[, t])) =", sum(is.na(wt[, t])), ", sum(wt[, t] =", sum(wt[, t]),"\n")
    W[, t] <- wt[, t] / sum(wt[,t])
    # resample particles by sampling parent particles according to weights
    A[, t] <- sample(1:npart, prob=W[1:npart, t], replace=T)

    # Resample particles for corresponding variables
    latent_var[, t,] <- latent_var[A[, t], t,]
    beta[t,] <- beta[t, A[, t]] #- needed for random walk on beta
    # beta_sd <- 5 * sd(beta[t,])
  } # end particle loop

  # Marginal likelihoods
  lik_values <- rep(NA, tend)
  for (t in 1:tend) {
    lik_values[t] <- log(sum(wt[1:npart, t])) # log-likelihoods
  }# averaged log likelihoods log(L/(npart^tend))
  loglik <- - tend * log(npart) + sum(lik_values)

  traj <- replicate(nstatevar + 1, rep(NA, tend), simplify=F)
  names(traj) <- c(names(y), "beta")
  # Latent variables sampled from the filtered distribution (backward sampling)
  if (backward_sampling) {
    loc <- rep(NA, tend)
    # loc_tend <-  sample(1:npart, prob = W[1:npart, tend], replace = T)
    loc[tend] <- sample.int(npart, size=1, prob=W[, tend], replace=T)
    # particle for the last time step
    traj[["beta"]][tend] <- beta[tend, loc[tend]]
    for (nm in names(y)) {
      traj[[nm]][tend] <- latent_var[loc[tend], tend, nm]
    }
    # lapply(names(y), function(x) traj[[x]][tend] <- latent_var[loc[tend], tend, x])
    # update backward
    for (i in seq(tend, 2, -1)) {
      loc[i-1] <- A[loc[i], i]
      traj[["beta"]][i-1] <- beta[i-1, loc[i-1]]
      for (nm in names(y)) {
        traj[[nm]][i-1] <- latent_var[loc[i-1], i-1, nm]
      }
    }
  }

  return (list(trace = traj, lik_marginal = lik_values,
               lik_overall_average = loglik,
               latent_var_filtered = latent_var,
               beta_filtered = beta,
               W = W, A = A))
}
