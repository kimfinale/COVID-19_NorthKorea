get_sum_stat_sepiar2 <- function(pars) {
  # library(Rcpp)
  # sourceCpp("src/sepiar_stoch.cpp")
  params <- params_init()
  # pars[1] gives the day (how many days earlier) when introduction of three cases
  p1 <- pars[1]
  params$ndays <- p1 + params$obslength # number of days for output
  params$R0 <- pars[2] # number of days for output
  params$day_intervention <- p1 + pars[3] # number of days for output
  params$R0_2 <- pars[4] # number of days for output

  out <- sepiar_stoch(params)
  day_filter <- seq(1, by = round(1/params$tau), length.out = params$ndays)
  out <- out[day_filter, "CI"]
  daily_inc <- diff(out)

  res <- c(sum(daily_inc[1:p1]), daily_inc[(p1+1):(p1+params$obslength)])

  return(res)
}
