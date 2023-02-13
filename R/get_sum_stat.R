get_sum_stat <- function(pars) {
  daily_inc <- daily_incidence(pars=pars)
  # res <- daily_incidence(pars=pars)
  # daily_inc <- res[, "CI"]
  p1 <- pars[1]
  res <- c(sum(daily_inc[1:p1]), daily_inc[(p1+1):length(daily_inc)])

  return(res)
}
