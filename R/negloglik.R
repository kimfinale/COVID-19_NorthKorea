negloglik <- function(pars, obs, dist="pois", size=50, ...) {

  inc <- incidence(pars=pars[1:4], ...)

  if (ncol(inc) > 1) {
    stop("More than one column in the daily incidence output")
  }
  p1 <- round(pars[1])
  model <- c(sum(inc[1:p1, 1]), inc[(p1+1):nrow(inc), 1])
  # model <- c(sum(inc[1:(p1+1), 1]), inc[(p1+2):nrow(inc), 1])
  # model <- inc[, 1]
  if (length(pars) >= 5) {
    size = pars[5] # NOTE: I ASSUME that the size parameter comes from the fifth position
  }
  # if the length of obs is different that of model, set the two to have the same size
  # this is useful when we use a fraction of observations to fit the model
  model <- model[1:length(obs)]
  if (any(model < 0)) {
    ll <- -Inf
  } else {
    if (dist == "pois") {
      ll <- sum(dpois(obs, model, log=TRUE), na.rm=TRUE)
    } else if(dist == "negbin") { # negative binomial distribution
      ll <- sum(dnbinom(obs, mu=model, size=size, log=TRUE), na.rm=T)
    }
    else{

    }
  }
  return(-ll)
}
