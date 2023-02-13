params_init <- function(erlang=FALSE){
  I0 <- 3 # since SEPIAR model is a stochastic model stochastic die-out possible
  y0 <- c(S=25970000, E=0, P=0, A=0, I=I0, R=0, CE=0, CI=0)
  y0 <- round(y0) # make sure that the y0 are integers

  params <- list() # input parameters for the SEIR model
  # epsilon and gamma from Kim et al. (2021)
  params$epsilon <- 1/3 # mean latent period = 1/epsilon
  params$delta <- 1/5.2 # mean incubation period = 1/delta
  params$gamma <- 1/6.5 # mean infectious period = 1/gamma

  params$R0 <- 3.0
  params$fA <- 0.306 # fraction of asymptomatic state
  params$bP <- 1 # relative infectiousness of pre-symptomatic state
  params$bA <- 1 # relative infectiousness of asymptomatic state
  params$tau <- 1 # time step size
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
    y0 <- c(S=25970000, E1=0, E2=0, P1=0, P2=0, A1=0, A2=0, I1=I0, I2=0, R=0, CE=0, CI=0)
    y0 <- round(y0) # make sure that the y0 are integers
    params$init <- NULL
    params$init$S <- y0[["S"]]
    params$init$E1 <- y0[["E1"]]
    params$init$E2 <- y0[["E2"]]
    params$init$P1 <- y0[["P1"]]
    params$init$P2 <- y0[["P2"]]
    params$init$I1 <- y0[["I1"]]
    params$init$I2 <- y0[["I2"]]
    params$init$A1 <- y0[["A1"]]
    params$init$A2 <- y0[["A2"]]
    params$init$R <- y0[["R"]]
    params$init$CE <- y0[["CE"]]
    params$init$CI <- y0[["CI"]]
  }

  return(params)
}
