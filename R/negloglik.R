negloglik <- function(pars, obs, dist="pois", size=50, ...) {

  inc <- incidence(pars=pars, ...)

  if (ncol(inc) > 1) {
    stop("More than one column in the daily incidence output")
  }
  p1 <- round(pars[1])
  model <- c(sum(inc[1:p1, 1]), inc[(p1+1):nrow(inc), 1])
  # model <- c(sum(inc[1:(p1+1), 1]), inc[(p1+2):nrow(inc), 1])
  # model <- inc[, 1]
  if (any(model < 0)) {
    ll <- -Inf
  } else {
    if (dist == "pois") {
      ll <- sum(dpois(obs, model, log=TRUE), na.rm=TRUE)
    } else if(dist == "negbin") { # negative binomial distribution
      ll <- sum(dnbinom(obs, size=size, mu=model, log=TRUE), na.rm=T)
    }
    else{

    }
  }
  return(-ll)
}
