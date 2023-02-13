# params2 indicate the paramters that were estimated
run_sepiar <- function(params, model, params2) {
  # params <- params_init(erlang=erlang)
  # if (erlang) params$tau <- 0.2 # more compartments -> higher rates -> smaller time steps
  maxd <- max(round(params2[,1])) # maximum value
  incmat <- matrix(NA, nrow=(maxd + params$obslength), ncol=nrow(params2)) # to store symptomatic
  totincmat <- incmat # to store total infection

  for (i in 1:nrow(params2)) {
    # cat("i = ", i, "\n")
    eps <- maxd - round(params2[i, 1])
    params$ndays <- round(params2[i, 1]) + params$obslength + 1
    params$R0 <- params2[i, 2]
    params$day_intervention <- round(params2[i, 1]) + round(params2[i, 3]) # number of days for output
    params$R0_2 <- params2[i, 4]

    out <- model(params)

    day_filter <- seq(1, by=ceiling(1/params$tau), length.out=params$ndays)
    incmat[(eps+1):nrow(incmat), i] <- diff(out[day_filter, "CI"]) # daily symptomatic
    totincmat[(eps+1):nrow(incmat), i] <- diff(out[day_filter, "CE"])
  }

  inc_summary <-
    as.data.frame(t(apply(incmat, 1, quantile,
                          probs=c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T)))
  totinc_summary <-
    as.data.frame(t(apply(totincmat, 1, quantile,
                          probs=c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T)))

  inc_summary$day <- 0:(nrow(inc_summary)-1)
  inc_summary$date <-
    seq(as.Date("2022-05-13") - maxd,
        as.Date("2022-05-13") + params$obslength - 1, by="day")

  totinc_summary$day <- inc_summary$day
  totinc_summary$date <- inc_summary$date

  return(list(infection=totinc_summary, symptomatic=inc_summary))
}

