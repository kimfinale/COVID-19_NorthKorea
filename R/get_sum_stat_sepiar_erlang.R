get_sum_stat_sepiar_erlang <- function(pars) {
  params <- params_init(erlang=TRUE)
  params$tau <- 0.2 #
  # pars[1] gives the day (how many days earlier) when introduction of three cases
  p1 <- round(pars[1]) # integer day makes it clear to model
  params$ndays <- p1 + params$obslength + 1 # number of days for output
  params$R0 <- pars[2] # number of days for output
  # previous modeling analyses revealed that almost 0 day (May 13), the intervention started
  # therefore, arbitrary negative number, -10, were added to explore
  # the possibility of intervention started earlier
  params$day_intervention <- p1 + pars[3] # number of days for output
  params$R0_2 <- pars[4] # number of days for output

  out <- sepiar_erlang_stoch(params)
  day_filter <- seq(1, by = round(1/params$tau), length.out = params$ndays)
  out <- out[day_filter, "CI"]
  daily_inc <- diff(out)

  res <- c(sum(daily_inc[1:p1]), daily_inc[(p1+1):(p1+params$obslength)])

  return(res)
}
