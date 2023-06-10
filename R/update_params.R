update_params <- function(pars=NULL) {
  params = PARAMETERS # PARAMETERS is a global variable holding baseline parameters
  # this is not the most efficient way
  params[["Day1"]] = pars[1]
  params[["R0"]] = pars[2]
  params[["Day2"]] = pars[3]
  params[["R0_int"]] = pars[4]

  # simulation times changes by Day 1 (introduction of the index case)
  params[["ndays"]] <- round(params[["Day1"]]) + params[["obs_length"]] + 1
  # day of intervention is counted by day
  params[["day_intervention"]] <- round(params[["Day1"]]) + round(params[["Day2"]])

  return (params)
}
