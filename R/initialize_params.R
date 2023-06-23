initialize_params <- function(...) {
  params <- list() # input parameters for the model
  params$population <- 25970000
  params$susceptible <- params$population - 1
  params$exposed <- 0
  params$asymptomatic <- 0
  params$presymptomatic <- 0
  params$symptomatic <- 1
  params$recovered <- 0
  params$dead <- 0
  params$cumul_infected <- 0
  params$cumul_symptomatic <- 0


  params$Day1 <- 0 # introduction happened "Day1" days before May 12 (reported)
  params$Day2 <- 4 # intervention "Day2" days after May 12
# 3 parameters to be estimated
  params$prop_inf <- 0.01
  params$R0 <- 6.0
  params$R0_int <- 0.6 # R0 after intervention is in place

  params$model <- seapird_euler
  # params$erlang <- FALSE
  # params$region <- "overall"
  params$epsilon <- 1/2.7 # mean latent period = 1/epsilon. Jiang (2023) Chinese doi:10.3760/cma.j.cn112150-20220926-00925
  params$delta <- 1/3.4 # mean incubation period = 1/delta. Wu (2022) JAMA Network Open doi:10.1001/jamanetworkopen.2022.28008
  params$gamma <- 1/5 # mean infectious period = 1/gamma. Takahashi (2022) EID doi:10.3201/eid2805.220197
  # this is solely based on a reasonable assumption that
  # case fatality ratio was measured as deaths observed within 30 days and
  # Gamma distribution with shape = 2 and rate = 2 * eta
  # Strasser (2022) JAMA Network Open doi:10.1001/jamanetworkopen.2022.38354
  params$eta <- 1 / 14 # 1 / mean delay from symptom onset to death

  params$fA <- 0.255 # fraction of asymptomatic state Yu (2022) JMW doi:10.1002/jmv.28066
  params$bP <- 1 # relative infectiousness of pre-symptomatic state
  params$bA <- 1 # relative infectiousness of asymptomatic state
  params$cfr <- 0.0304 # case fatality ratio Wang(2023) JMV doi:10.1002/jmv.28118

  params$tau <- 0.1 # time step size
  params$ndays <- 100.0 # number of days for output
  params$day_intervention <- 100.0

  # North Korea data have >80 observations, ignoring some zero incidence at the end
  params$obs_length <- 74


  # update parameters
  par = list(...)

  nms = names(par)
  for (nm in nms) {
    params[[nm]] = par[[nm]]
  }

  fa <- params[["fA"]]
  params[["asymptomatic"]] <- params[["symptomatic"]] * (fa)/(1-fa)
  infecteds <- params[["asymptomatic"]] + params[["symptomatic"]]
  params[["susceptible"]] <- params[["population"]] - infecteds
  # simulation times changes by Day 1 (introduction of the virus)
  # + 1 because observation is assumed to be the difference between the time steps
  params[["ndays"]] <- round(params[["Day1"]]) + params[["obs_length"]] + 1
  # day of intervention is counted by days from introduction
  params[["day_intervention"]] <- round(params[["Day1"]]) + round(params[["Day2"]])

  return (params)
}
