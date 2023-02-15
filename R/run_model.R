run_model <- function(pars, var) {
  max_intro <- max(round(pars[,1])) # value for the earliest introduction

  reslist <- vector('list', length(var))
  names(reslist) <- var

  # PARAMETERS$measure_var <- var

  nobs <- max_intro + PARAMETERS$obslength
  np <- nrow(pars)

  for(j in 1:length(var)) {
    reslist[[j]] <- data.frame(matrix(NA, nrow=nobs, ncol=np))
  }

  for (i in 1:np) {
    res <- daily_incidence(pars = pars[i,])
    for(j in 1:length(var)) {
      reslist[[j]][(nobs-nrow(res)+1):nobs, i] <- res[, var[j]]
    }
  }

  return(reslist)
}
