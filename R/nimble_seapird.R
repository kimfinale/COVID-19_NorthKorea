nimble_seapird <- nimbleFunction(
  run = function(R0=double(0),
                 R_int=double(0),
                 pop=double(0),
                 inf0=double(0, default=1.0),
                 day1=double(0, default=0),
                 day2=double(0),
                 full_inc=integer(0,default=0)) {

  Day1 = day1 # introduction happened "Day1" days before May 12 (reported)
  Day2 = day2 # intervention "Day2" days after May 12
  epsilon = 1/2.7 # mean latent period = 1/epsilon. Jiang (2023) Chinese doi:10.3760/cma.j.cn112150-20220926-00925
  delta = 1/3.4 # mean incubation period = 1/delta. Wu (2022) JAMA Network Open doi:10.1001/jamanetworkopen.2022.28008
  gamma = 1/5 # mean infectious period = 1/gamma. Takahashi (2022) EID doi:10.3201/eid2805.220197
  eta = 1/14 # 1 / mean delay from symptom onset to death
  fA = 0.255 # fraction of asymptomatic state Yu (2022) JMW doi:10.1002/jmv.28066
  bP = 1     # relative infectiousness of pre-symptomatic state
  bA = 1     # relative infectiousness of asymptomatic state
  cfr = 0.0304 # case fatality ratio Wang(2023) JMV doi:10.1002/jmv.28118
  obs_length = 74                  # May 14 - July 26
  R0 = R0 # basic reproduction number
  R_int = R_int # R0 when intervention is in place
  rate_P_I = 1 / (1/delta - 1/epsilon) # transition rate from P to I
  # // 'dur_infect' is a temporary variable reflecting the fraction (fA) of going
  # // to A rather than P (which becomes I), the duration of P, A, and I
  # // stages, and the infectivity of P (ie, bP) and A (ie, bA) relative to I.
  dur_infect = fA * bA * (1/rate_P_I + 1/gamma) + (1 - fA) * (bP/rate_P_I + 1/gamma)
  beta = R0 / dur_infect # // 'dur_infect' multiplied by beta gives R0
  ndays = floor(Day1) + obs_length + 1
  # day of intervention is counted by days from introduction
  day_intervention = floor(Day1) + Day2
  # // Calculate the number of events for each step, update state vectors
  S = rep(0, ndays)
  E = rep(0, ndays)
  A = rep(0, ndays)
  P = rep(0, ndays)
  I = rep(0, ndays)
  R = rep(0, ndays)
  D = rep(0, ndays)
  CI = rep(0, ndays)

  # initial conditions
  It = inf0                         # symptomatic
  At = It * (fA)/(1-fA)             # asymptomatic
  infecteds = It + At               # total infecteds
  St = pop - infecteds              # susceptible
  Et = 0                            # exposed
  Pt = 0                            # pre-symptomatic
  Rt = 0                            # recovered
  Dt = 0                            # deaths
  CIt = 0                           # cumulative symptomatic

  tau <- 0.1 # time step size in terms of a day
## -----------------------------------------------------------------
  # time is treated as real numbers and
  # digits are used to create steps before day 1
  #
  time_before_day1 = Day1 - floor(Day1)
  steps_before_day1 = ceiling(time_before_day1 / tau)
  for(i in 1:steps_before_day1) {# // Essentially Euler method was implemented
    #// State Equations
    Nt = St + Et + Pt + It + At + Rt
    foi = beta * (bP*Pt + bA*At + It) / Nt # force of infection
    new_infections = St * foi * tau
    EtoP = Et * epsilon * (1-fA) * tau
    EtoA = Et * epsilon * fA * tau
    new_symptoms = Pt * rate_P_I * tau
    ItoR = It * (1-cfr) * gamma * tau
    AtoR = At * 1 / (1/gamma + 1/rate_P_I) * tau
    new_deaths = It * cfr * eta * tau
    # Calculate the change in each state variable
    dS = - new_infections
    dE = new_infections - EtoP - EtoA
    dP = EtoP - new_symptoms
    dI = new_symptoms - ItoR - new_deaths
    dA = EtoA - AtoR
    dR = ItoR + AtoR
    dD = new_deaths
    # // Update next timestep
    St = St + dS
    Et = Et + dE
    At = At + dA
    Pt = Pt + dP
    It = It + dI
    Rt = Rt + dR
    Dt = Dt + dD
    CIt = CIt + new_symptoms # cumulative symptomatic
  }
## -----------------------------------------------------------------
  S[1] = St
  E[1] = Et
  A[1] = At
  P[1] = Pt
  I[1] = It
  R[1] = Rt
  D[1] = Dt
  CI[1] = CIt

  stepsperday = ceiling(1/tau)
  telapsed = 0

  for (day in 2:ndays) {
    for(step in 1:stepsperday) {# // Essentially Euler method was implemented
      telapsed = telapsed + tau # time tracking
      if (telapsed >= day_intervention) {
        beta = R_int / dur_infect;
      }
      #// State Equations
      Nt = St + Et + Pt + It + At + Rt
      foi = beta * (bP*Pt + bA*At + It) / Nt # force of infection
      new_infections = St * foi * tau
      EtoP = Et * epsilon * (1-fA) * tau
      EtoA = Et * epsilon * fA * tau
      new_symptoms = Pt * rate_P_I * tau
      ItoR = It * (1-cfr) * gamma * tau
      AtoR = At * 1 / (1/gamma + 1/rate_P_I) * tau
      new_deaths = It * cfr * eta * tau

      # Calculate the change in each state variable
      dS = - new_infections
      dE = new_infections - EtoP - EtoA
      dP = EtoP - new_symptoms
      dI = new_symptoms - ItoR - new_deaths
      dA = EtoA - AtoR
      dR = ItoR + AtoR
      dD = new_deaths
      # // Update next timestep
      St = St + dS
      Et = Et + dE
      At = At + dA
      Pt = Pt + dP
      It = It + dI
      Rt = Rt + dR
      Dt = Dt + dD
      CIt = CIt + new_symptoms # cumulative symptomatic
    }

    S[day] = St
    E[day] = Et
    A[day] = At
    P[day] = Pt
    I[day] = It
    R[day] = Rt
    D[day] = Dt
    CI[day] = CIt
  }


  inc = CI[(2+floor(day1)):ndays] - CI[(1+floor(day1)):(ndays-1)]

  if (full_inc > 1e-6) {
    inc = CI[2:ndays] - CI[1:(ndays-1)]
  }
  return(inc)
  returnType(double(1))
}
)
