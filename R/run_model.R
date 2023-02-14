run_model <- function(pars, outvar) {
  max_intro <- max(round(pars[,1])) # value for the earliest introduction
  reslist <- vector('list', length(outvar))
  PARAMETERS$measure_var <- outvar

  nobs <- max_intro + PARAMETERS$obslength
  np <- nrow(pars)

  for(i in 1:length(outvar)) {
    reslist[[i]] <- data.frame(matrix(NA, nrow=nobs, ncol=np))
  }
  names(reslist) <- outvar
  for (i in 1:np) {
    res <- daily_incidence(pars = pars[i,])
    for(j in 1:length(outvar)) {
      reslist[[j]][(nobs-nrow(res)+1):nobs, i] <- res[, outvar[j]]
    }
  }

  return(reslist)
}
