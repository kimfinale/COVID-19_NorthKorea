negloglik <- function(pars, dat, error_dist="pois") {
  daily_inc <- daily_incidence(pars=pars)
  if (ncol(daily_inc) > 1){
    stop("More than one column in the daily incidence output")
  }
  p1 <- round(pars[1])
  model <- c(sum(daily_inc[1:p1, 1]),
             daily_inc[(p1+1):nrow(daily_inc), 1])
  if (any(model < 0)) {
    ll <- -Inf
  } else {
    if (error_dist == "pois") {
      ll <- sum(dpois(dat, model, log=TRUE))
    } else { # negative binomial distribution
      ll <- sum(dnbinom(dat, size=pars[5], mu=model, log=TRUE))
    }
  }
  return(-ll)
}
