#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector reulermultinom(double size, NumericVector rate, double dt) {
  int ncol = rate.size();
  NumericVector trans(ncol); // transition events
  double p = sum(rate); //total event rate
  double tmpp = p;
  double tmpsize = R::rbinom(size, (1 - exp(-tmpp*dt))); // total number of events
  for (int k = 0; k < (ncol-1); k++) {
    double tr = R::rbinom(tmpsize, rate(k)/p);
    trans(k) = tr;
    tmpsize -= tr;
    tmpp -= rate(k);
  }
  trans(ncol-1) = tmpsize;
  return(trans);
}

// the model framework adopted from https://gallery.rcpp.org/articles/epidemiological-compartment-model/
// [[Rcpp::export]]
List sepiar_stoch(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  NumericVector S(nsteps);
  NumericVector E(nsteps);
  NumericVector P(nsteps);
  NumericVector I(nsteps);
  NumericVector A(nsteps); //asymptomatic
  NumericVector R(nsteps); // recovered

  NumericVector CE(nsteps); // cumulative infection
  NumericVector CI(nsteps); // cumulative symptomatic

  NumericVector time(nsteps);// fill time w/zeros

  // initial values
  List init = params["init"];
  S(0) = init["S"];
  E(0) = init["E"];
  P(0) = init["P"];
  I(0) = init["I"];
  A(0) = init["A"]; //asymptomatic
  R(0) = init["R"]; // recovered
  CE(0) = init["CE"]; // cumulative infection
  CI(0) = init["CI"]; // cumulative symptomatic

  double epsilon = params["epsilon"]; // 1 / latent period
  double delta = params["delta"]; // 1 / incubation period
  double gamma = params["gamma"]; // 1 / recovery period

  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double rate_P_I = 1 / (1/delta - 1/epsilon);
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A compared to I
  double bP = params["bP"]; // relative infectivity of P compared to I
  double R0 = params["R0"];
  double R0_2 = params["R0_2"]; // fraction of R0 when intervention is in place
  double day_intervention = params["day_intervention"];

  double beta = R0 / (bP/rate_P_I + bA*fA/gamma + (1-fA)/gamma);

  // Rprintf("the value of delta : %.2f \n", delta);
  // Rprintf("the value of nsteps : %.1d \n", nsteps);
  // Rprintf("the value of delta : %.2f \n", delta);
  // Rprintf("the value of epsilon : %.2f \n", epsilon);
  // Rprintf("the value of gamma : %.2f \n", gamma);
  // Rprintf("the value of rate_P_I : %.2f \n", rate_P_I);
  //
  NumericVector P_rates = {rate_P_I*(1-fA), rate_P_I*fA};

  // Rprintf("the value of P_rates[0] : %.2f \n", P_rates[0]);
  // Rprintf("the value of P_rates[1] : %.2f \n", P_rates[1]);

  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < nsteps - 1; istep++) {

    if (istep*tau >= day_intervention) {
      beta = R0_2 / (bP/rate_P_I + bA*fA/gamma + (1-fA)/gamma);

    }

    double iS = S[istep];
    double iE = E[istep];
    double iP = P[istep];
    double iI = I[istep];
    double iA = A[istep];
    double iR = R[istep];
    double iCE = CE[istep];
    double iCI = CI[istep];

    // State Equations
    double N = iS + iE + iP + iI + iA + iR;
    double inf = bP*iP + bA*iA + iI;
    double foi = beta * inf / N;

    // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
    double new_infection = R::rbinom(iS, 1 - exp(-foi * tau));
    double EtoP = R::rbinom(iE, 1 - exp(- epsilon * tau)); //

    // Rprintf("the value of new_infection : %.1f \n", new_infection);
    // Rprintf("the value of EtoP : %.1f \n", EtoP);

    NumericVector from_P = reulermultinom(iP, P_rates, tau);
    double PtoI = from_P(0); //
    double PtoA = from_P(1); //

    // Rprintf("the value of PtoI : %.1f \n", PtoI);
    // Rprintf("the value of PtoA : %.1f \n", PtoA);

    double ItoR = R::rbinom(iI, 1 - exp(- gamma * tau));
    double AtoR = R::rbinom(iA, 1 - exp(- gamma * tau)); //

    // Rprintf("the value of ItoR : %.1f \n", ItoR);
    // Rprintf("the value of AtoR : %.1f \n", AtoR);

    // Calculate the change in each state variable
    double dS = - new_infection;
    double dE = new_infection - EtoP;
    double dP = EtoP - PtoI - PtoA;
    double dI = PtoI - ItoR;
    double dA = PtoA - AtoR;
    double dR = ItoR + AtoR;

    // Update next timestep
    S[istep + 1] = iS + dS;
    E[istep + 1] = iE + dE;
    P[istep + 1] = iP + dP;
    I[istep + 1] = iI + dI;
    A[istep + 1] = iA + dA;
    R[istep + 1] = iR + dR;
    CE[istep + 1] = iCE + new_infection;// cumulative infection
    CI[istep + 1] = iCI + PtoI;// cumulative symptomatic
    time[istep + 1] = (istep + 1) * tau;// time in fractional years

  }
  // Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S") = S,
    Named("E") = E,
    Named("P") = P,
    Named("I") = I,
    Named("A") = A,
    Named("R") = R,
    Named("CE") = CE,
    Named("CI") = CI);

  return sim;
}


// the model framework adopted from https://gallery.rcpp.org/articles/epidemiological-compartment-model/
// [[Rcpp::export]]
List sepiar_erlang_stoch(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  NumericVector S(nsteps);
  NumericVector E1(nsteps);
  NumericVector E2(nsteps);
  NumericVector P1(nsteps);
  NumericVector P2(nsteps);
  NumericVector I1(nsteps);
  NumericVector I2(nsteps);
  NumericVector A1(nsteps); //asymptomatic
  NumericVector A2(nsteps);
  NumericVector R(nsteps); // recovered

  NumericVector CE(nsteps); // cumulative infection
  NumericVector CI(nsteps); // cumulative symptomatic

  NumericVector time(nsteps);// fill time w/zeros

  // initial values
  List init = params["init"];
  S(0) = init["S"];
  E1(0) = init["E1"];
  E2(0) = init["E2"];
  P1(0) = init["P1"];
  P2(0) = init["P2"];
  I1(0) = init["I1"];
  I2(0) = init["I2"];
  A1(0) = init["A1"]; //asymptomatic
  A2(0) = init["A2"]; //asymptomatic
  R(0) = init["R"]; // recovered
  CE(0) = init["CE"]; // cumulative infection
  CI(0) = init["CI"]; // cumulative symptomatic

  double epsilon = params["epsilon"]; // 1 / latent period
  double delta = params["delta"]; // 1 / incubation period
  double gamma = params["gamma"]; // 1 / recovery period

  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double rate_P_I = 1 / (1/delta - 1/epsilon);
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A compared to I
  double bP = params["bP"]; // relative infectivity of P compared to I
  double R0 = params["R0"];
  double R0_2 = params["R0_2"]; // fraction of R0 when intervention is in place
  double day_intervention = params["day_intervention"];

  double beta = R0 / (bP/rate_P_I + bA*fA/gamma + (1-fA)/gamma);

  NumericVector P2_rates = {rate_P_I*2*(1-fA), rate_P_I*2*fA};

  // Rprintf("the value of P_rates[0] : %.2f \n", P_rates[0]);
  // Rprintf("the value of P_rates[1] : %.2f \n", P_rates[1]);

  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < nsteps - 1; istep++) {

    if (istep*tau >= day_intervention) {
      beta = R0_2 / (bP/rate_P_I + bA*fA/gamma + (1-fA)/gamma);

    }

    double iS = S[istep];
    double iE1 = E1[istep];
    double iE2 = E1[istep];
    double iP1 = P1[istep];
    double iP2 = P2[istep];
    double iI1 = I1[istep];
    double iI2 = I2[istep];
    double iA1 = A1[istep];
    double iA2 = A2[istep];
    double iR = R[istep];
    double iCE = CE[istep];
    double iCI = CI[istep];

    // State Equations
    double N = iS + iE1 + iE2 + iP1 + iP2 + iI1 + iI2 + iA1 + iA2 + iR;
    double inf = bP*(iP1 + iP2) + bA*(iA1 + iA2) + iI1 + iI2;
    double foi = beta * inf / N;

    // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
    double new_infection = R::rbinom(iS, 1 - exp(-foi * tau));
    double E1toE2 = R::rbinom(iE1, 1 - exp(- epsilon * 2 * tau)); //
    double E2toP1 = R::rbinom(iE2, 1 - exp(- epsilon * 2 * tau)); //
    double P1toP2 = R::rbinom(iP1, 1 - exp(- rate_P_I * 2 * tau)); //

    NumericVector from_P2 = reulermultinom(iP2, P2_rates, tau);
    double P2toI1 = from_P2(0); //
    double P2toA1 = from_P2(1); //

    double A1toA2 = R::rbinom(iA1, 1 - exp(- gamma * 2 * tau)); //
    double A2toR = R::rbinom(iA2, 1 - exp(- gamma * 2 * tau)); //
    double I1toI2 = R::rbinom(iI1, 1 - exp(- gamma * 2 * tau));
    double I2toR = R::rbinom(iI2, 1 - exp(- gamma * 2 * tau));


    // Calculate the change in each state variable
    double dS = - new_infection;
    double dE1 = new_infection - E1toE2;
    double dE2 = E1toE2 - E2toP1;
    double dP1 = E2toP1 - P1toP2;
    double dP2 = P1toP2 - P2toI1 - P2toA1;
    double dI1 = P2toI1 - I1toI2;
    double dI2 = I1toI2 - I2toR;
    double dA1 = P2toA1 - A1toA2;
    double dA2 = A1toA2 - A2toR;
    double dR = I2toR + A2toR;

    // Update next timestep
    S[istep + 1] = iS + dS;
    E1[istep + 1] = iE1 + dE1;
    E2[istep + 1] = iE2 + dE2;
    P1[istep + 1] = iP1 + dP1;
    P2[istep + 1] = iP2 + dP2;
    I1[istep + 1] = iI1 + dI1;
    I2[istep + 1] = iI2 + dI2;
    A1[istep + 1] = iA1 + dA1;
    A2[istep + 1] = iA2 + dA2;
    R[istep + 1] = iR + dR;
    CE[istep + 1] = iCE + new_infection;// cumulative infection
    CI[istep + 1] = iCI + P2toI1;// cumulative symptomatic
    time[istep + 1] = (istep + 1) * tau;// time in fractional years

  }
  // Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S") = S,
    Named("E1") = E1,
    Named("E2") = E2,
    Named("P1") = P1,
    Named("P2") = P2,
    Named("I1") = I1,
    Named("I2") = I2,
    Named("A1") = A1,
    Named("A2") = A2,
    Named("R") = R,
    Named("CE") = CE,
    Named("CI") = CI);

  return sim;
}



// the model framework adopted from https://gallery.rcpp.org/articles/epidemiological-compartment-model/
// [[Rcpp::export]]
List sepiar_erlang_euler(List params) {
  double tau = params["tau"]; // time step size
  double ndays = params["ndays"]; // number of days for output
  int nsteps = ceil(ndays / tau) + 1;
  // vectors for state variables
  NumericVector S(nsteps);
  NumericVector E1(nsteps);
  NumericVector E2(nsteps);
  NumericVector P1(nsteps);
  NumericVector P2(nsteps);
  NumericVector I1(nsteps);
  NumericVector I2(nsteps);
  NumericVector A1(nsteps); //asymptomatic
  NumericVector A2(nsteps);
  NumericVector R(nsteps); // recovered

  NumericVector CE(nsteps); // cumulative infection
  NumericVector CI(nsteps); // cumulative symptomatic

  NumericVector time(nsteps);// fill time w/zeros

  // initial values
  List init = params["init"];
  S(0) = init["S"];
  E1(0) = init["E1"];
  E2(0) = init["E2"];
  P1(0) = init["P1"];
  P2(0) = init["P2"];
  I1(0) = init["I1"];
  I2(0) = init["I2"];
  A1(0) = init["A1"]; //asymptomatic
  A2(0) = init["A2"]; //asymptomatic
  R(0) = init["R"]; // recovered
  CE(0) = init["CE"]; // cumulative infection
  CI(0) = init["CI"]; // cumulative symptomatic

  double epsilon = params["epsilon"]; // 1 / latent period
  double delta = params["delta"]; // 1 / incubation period
  double gamma = params["gamma"]; // 1 / recovery period

  // 1 / rate from pre-symptomatic (P) state to infectious (I) states
  double rate_P_I = 1 / (1/delta - 1/epsilon);
  double fA = params["fA"]; // fraction asymptomatic
  double bA = params["bA"]; // relative infectivity of A compared to I
  double bP = params["bP"]; // relative infectivity of P compared to I
  double R0 = params["R0"];
  double R0_2 = params["R0_2"]; // fraction of R0 when intervention is in place
  double day_intervention = params["day_intervention"];

  double beta = R0 / (bP/rate_P_I + bA*fA/gamma + (1-fA)/gamma);

  // Rprintf("the value of P_rates[0] : %.2f \n", P_rates[0]);
  // Rprintf("the value of P_rates[1] : %.2f \n", P_rates[1]);

  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < nsteps - 1; istep++) {

    if (istep*tau >= day_intervention) {
      beta = R0_2 / (bP/rate_P_I + bA*fA/gamma + (1-fA)/gamma);

    }

    double iS = S[istep];
    double iE1 = E1[istep];
    double iE2 = E1[istep];
    double iP1 = P1[istep];
    double iP2 = P2[istep];
    double iI1 = I1[istep];
    double iI2 = I2[istep];
    double iA1 = A1[istep];
    double iA2 = A2[istep];
    double iR = R[istep];
    double iCE = CE[istep];
    double iCI = CI[istep];

    // State Equations
    double N = iS + iE1 + iE2 + iP1 + iP2 + iI1 + iI2 + iA1 + iA2 + iR;
    double inf = bP*(iP1 + iP2) + bA*(iA1 + iA2) + iI1 + iI2;
    double foi = beta * inf / N;

    // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
    double new_infection = iS * foi * tau;
    double E1toE2 = iE1 * epsilon * 2 * tau; //
    double E2toP1 = iE2 * epsilon * 2 * tau; //
    double P1toP2 = iP1 * rate_P_I * 2 * tau; //

    double P2toI1 = iP2 * rate_P_I * (1-fA) * 2 * tau; //
    double P2toA1 = iP2 * rate_P_I * fA * 2 * tau; //

    double A1toA2 = iA1 * gamma * 2 * tau; //
    double A2toR = iA2 * gamma * 2 * tau; //
    double I1toI2 = iI1 * gamma * 2 * tau;
    double I2toR = iI2 * gamma * 2 * tau;


    // Calculate the change in each state variable
    double dS = - new_infection;
    double dE1 = new_infection - E1toE2;
    double dE2 = E1toE2 - E2toP1;
    double dP1 = E2toP1 - P1toP2;
    double dP2 = P1toP2 - P2toI1 - P2toA1;
    double dI1 = P2toI1 - I1toI2;
    double dI2 = I1toI2 - I2toR;
    double dA1 = P2toA1 - A1toA2;
    double dA2 = A1toA2 - A2toR;
    double dR = I2toR + A2toR;

    // Update next timestep
    S[istep + 1] = iS + dS;
    E1[istep + 1] = iE1 + dE1;
    E2[istep + 1] = iE2 + dE2;
    P1[istep + 1] = iP1 + dP1;
    P2[istep + 1] = iP2 + dP2;
    I1[istep + 1] = iI1 + dI1;
    I2[istep + 1] = iI2 + dI2;
    A1[istep + 1] = iA1 + dA1;
    A2[istep + 1] = iA2 + dA2;
    R[istep + 1] = iR + dR;
    CE[istep + 1] = iCE + new_infection;// cumulative infection
    CI[istep + 1] = iCI + P2toI1;// cumulative symptomatic
    time[istep + 1] = (istep + 1) * tau;// time in fractional years

  }
  // Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S") = S,
    Named("E1") = E1,
    Named("E2") = E2,
    Named("P1") = P1,
    Named("P2") = P2,
    Named("I1") = I1,
    Named("I2") = I2,
    Named("A1") = A1,
    Named("A2") = A2,
    Named("R") = R,
    Named("CE") = CE,
    Named("CI") = CI);

  return sim;
}
