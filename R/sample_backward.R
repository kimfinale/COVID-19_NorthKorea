# backward sampling
sample_backward <- function (y, 
                             latent_var, 
                             varname_est = "beta", 
                             tend, 
                             npart, 
                             weight){
  # one trajectory for each variable plus beta we estimate
  nstatevar <- length(y)
  traj <- replicate(nstatevar + 1, rep(NA, tend), simplify = F)
  names(traj) <- c(names(y), varname_est)
  # Latent variables sampled from the filtered distribution (backward sampling)
  loc <- rep(NA, tend)
  # loc_tend <-  sample(1:npart, prob = W[1:npart, tend], replace = T)
  loc[tend] <-  sample.int(npart, size = 1, prob = W[, tend], replace = T)
  
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
