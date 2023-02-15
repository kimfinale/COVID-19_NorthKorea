daily_incidence <- function(pars){
  # use the global variable PARAMETERS, which holds all the parameters
  # and initial conditions
  # pars only stores estimated parameters
  params <- PARAMETERS
  p1 <- round(pars[1]) # integer day makes it clear to model
  # +1 because the output comes as the daily difference
  params$ndays <- p1 + params$obslength # number of days for output
  params$R0 <- pars[2] # number of days for output
  # previous modeling analyses revealed that almost 0 day (May 13), the intervention started
  # therefore, arbitrary negative number, -10, were added to explore
  # the possibility of intervention started earlier
  params$day_intervention <- p1 + pars[3]
  params$R0_2 <- pars[4] # number of days for o

  out <- params$model(params)
  day_filter <- seq(1, by=round(1/params$tau), length.out=(params$ndays+1))
  mv <- params$measure_var
  out <- out[day_filter, mv, drop=FALSE]

  # if measured variables is more than one
  df <- data.frame(matrix(NA, nrow=(nrow(out)-1), ncol=length(mv)))
  names(df) <- mv
  for(i in 1:length(mv)){
    df[,i] <- diff(out[, mv[i]])
  }

  return(df)
}
