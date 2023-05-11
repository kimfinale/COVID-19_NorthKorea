initialize_params <- function(I0=3, pop=25970000, erlang=FALSE, region="overall"){
  # 25970000 North Korea population size
  # since SEPIAR model is a stochastic model stochastic die-out possible
  if(region != "overall"){
    pop <- get_population_size(region)
  }
  y0 <- c(S=pop-I0, E=0, P=0, A=0, I=I0, R=0, CE=0, CI=0)
  y0 <- round(y0) # make sure that the y0 are integers

  params <- list() # input parameters for the SEIR model
  params$measure_var <- "CI" # daily diff of CI gives daily symptomatic case
  params$model <- sepiar_stoch
  params$erlang <- erlang
  params$region <- region
  # epsilon and gamma from Kim et al. (2021)
  params$epsilon <- 1/3 # mean latent period = 1/epsilon
  params$delta <- 1/5.2 # mean incubation period = 1/delta
  params$gamma <- 1/6.5 # mean infectious period = 1/gamma

  params$R0 <- 3.0
  params$fA <- 0.306 # fraction of asymptomatic state
  params$bP <- 1 # relative infectiousness of pre-symptomatic state
  params$bA <- 1 # relative infectiousness of asymptomatic state
  params$tau <- 0.1 # time step size
  params$ndays <- 100.0 # number of days for output
  params$day_intervention <- 100.0
  params$R0_2 <- 1.0 # fraction of R
  # North Korea data have >80 observations, ignoring some zero incidence
  params$obslength <- 80

  params$init$S <- y0[["S"]]
  params$init$E <- y0[["E"]]
  params$init$P <- y0[["P"]]
  params$init$I <- y0[["I"]]
  params$init$A <- y0[["A"]]
  params$init$R <- y0[["R"]]
  params$init$CE <- y0[["CE"]]
  params$init$CI <- y0[["CI"]]

  if (erlang) {
    params$model <- sepiar_erlang_stoch
    params$tau <- 0.2 # time step size
    params$init <- NULL
    params$init$S <- pop-I0
    params$init$E1 <- 0
    params$init$E2 <- 0
    params$init$P1 <- 0
    params$init$P2 <- 0
    params$init$I1 <- I0
    params$init$I2 <- 0
    params$init$A1 <- 0
    params$init$A2 <- 0
    params$init$R <- 0
    params$init$CE <- 0
    params$init$CI <- 0
  }

  return(params)
}
