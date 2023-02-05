get_sum_stat_sepiar <- function(pars) {
  # library(Rcpp)
  # sourceCpp("src/sepiar_stoch.cpp")
  params <- params_init()
  # pars[1] gives the day (how many days earlier) when introduction of three cases
  p1 <- pars[1]
  params$ndays <- p1 + 5 # number of days for output
  params$R0 <- pars[2] # number of days for output
  
  # params$fR0 <- pars[3] # number of days for output
  # params$day_intervention <- p1 + pars[4] # number of days for output
  
  params$fR0 <- 1 # number of days for output
  params$day_intervention <- p1 + 5 # number of days for output
  
  out <- sepiar_stoch(params)
  day_filter <- seq(1, by = round(1/params$tau), length.out = params$ndays)
  out <- out[day_filter, "CI"]
  daily_inc <- diff(out)
  
  res <- c(sum(daily_inc[1:p1]), daily_inc[(p1+1):(p1+4)])
  
  # return(c(sum(daily_inc[1:p1]), sum(daily_inc[(p1+1):length(daily_inc)])))
  # return(c(sum(daily_inc[1:p1]), daily_inc[(p1+1):length(daily_inc)]))
  # res <- c(sum(daily_inc[1:p1]), daily_inc[(p1+1):(p1+78)]) / 
  # c(350000, 18000, 174440, 296180, 392920, 269510, 232880, 262270, 
  #   263370, 219030, 186090, 167650, 134510, 115970, 105500, 100460, 
  #   88520, 89500, 100710, 96020, 93180, 96610, 82160, 79100, 73780, 
  #   66680, 61730, 54610, 50860, 45540, 42810, 40060, 36710, 32810, 
  #   29910, 26010, 23160, 20360, 19310, 18820, 17250, 15260, 13100, 
  #   11010, 9610, 8920, 7300, 6710, 5980, 4730, 4570, 4100, 3540, 
  #   3030, 2500, 2140, 1950, 1630, 1590, 1460, 1240, 900, 770, 560, 
  #   500, 460, 430, 310, 250, 250, 170, 140, 120, 120, 50, 30, 18, 
  #   11, 3)
  # return(c(sum(daily_inc[1:p1]), daily_inc[(p1+1):length(daily_inc)]))
  return(res)
}
