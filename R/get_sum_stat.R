get_sum_stat <- function(pars) {
  daily_inc <- daily_incidence(pars=pars)
  if (ncol(daily_inc) > 1){
    stop("More than one column in the daily incidence output")
  }
  p1 <- round(pars[1])
  res <- c(sum(daily_inc[1:p1, 1]), daily_inc[(p1+1):nrow(daily_inc), 1])

  return(res)
}
