#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>

/*******************************************************************
 *
 * Algae is a simplified model for the growth of algae as described
 * by Weber (2012)
 *
 * The model provides no additional output on intermediary variables on
 * request
 *
 *******************************************************************/

/**
 * Allocate memory for global parameter array
 */
static double parms[16] = {0};
/**
 * Allocate memory for forcing function data
 *
 * Array values get updated by the ODE solver in every time step.
 */
static double forc[3] = {0};
/*
 * Constant helper value
 */
static double log_EC50 = 0;
/*
 * Define aliases
 */
// State variable aliases
#define A y[0]
#define Q y[1]
#define P y[2]
#define Dw y[3]

// Derivative aliases
#define dA ydot[0]
#define dQ ydot[1]
#define dP ydot[2]
#define dDw ydot[3]

// Parameter aliases
// Growth model parameters
#define mu_max parms[0]
#define m_max parms[1]
#define v_max parms[2]
#define k_s parms[3]
#define Q_min parms[4]
#define Q_max parms[5]
#define T_opt parms[6]
#define T_min parms[7]
#define T_max parms[8]
#define I_opt parms[9]
#define R_0 parms[14]
#define D parms[15]
// Toxicodynamic parameters
#define EC_50 parms[10]
#define b parms[11]
// Toxicokinetic parameters
#define kD parms[12]
#define dose_resp parms[13]  // 0 = logit, 1 = probit
// Forcings by environmental variables
#define Cw    forc[0]
#define I     forc[1]
#define T_act forc[2]

/**
 * Parameter initializer
 */
void algae_tktd_init(void (*odeparms)(int *, double *))
{
  int N = 16;
  odeparms(&N, parms);

  log_EC50 = log(EC_50);
}

/**
 * Forcings initializer
 */
void algae_tktd_forc(void (*odeforcs)(int *, double *))
{
  int N = 3;
  odeforcs(&N, forc);
}

/*
 * Functions for calculating reduction factors for growth rate
 */
// Temperature effect on mu_max
double f_temp_tktd(double temp)
{
  double f = (temp < T_opt) ? (temp - T_opt) / (T_min - T_opt) : (temp - T_opt) / (T_max - T_opt);
  return exp(-2.3 * pow(f, 2));
}

// Irradiance effect on mu_max
double f_I_tktd(double I_param)
{
  return I_param / I_opt * exp(1 - (I_param / I_opt));
}

// Phosphorus effect on mu_max
double f_Q_tktd(double Q_param, double A_param)
{
  return (1 - exp(-log(2)*((Q_param / (Q_min * A_param)) - 1)));
}

// Effect on internal phosphorus
double f_Q_P_tktd(double A_param, double Q_param, double P_param)
{
  return (((Q_max * A_param - Q_param) / ((Q_max - Q_min) * A_param)) * (P_param / (k_s + P_param)));
}

// Effect of concentration on growth using probit (OECD 54, page 67)
double f_Dw_probit(double Dw_param) {
  //return gsl_cdf_ugaussian_P(-b * (log(C_param) - log(EC_50)));
  //Function to calculate the cumulative distribution function (CDF) of the standard normal distribution
  double z = -b * (log(Dw_param) - log_EC50);
  return 0.5 * erfc(-z / sqrt(2));
}

// Effect of concentration on growth using logit (OECD 54, page 67)
double f_Dw_logit(double Dw_param)
{
  return 1 / (1 + exp(b * (log(Dw_param) - log_EC50)));
}

/**
 * Derivatives
 */
void algae_tktd_func(int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{

  // Biomass
  if(dose_resp == 0) {
    dA = (mu_max * f_temp_tktd(T_act) * f_I_tktd(I) * f_Q_tktd(Q, A) * f_Dw_logit(Dw) - m_max - D) * A;
  } else {
    dA = (mu_max * f_temp_tktd(T_act) * f_I_tktd(I) * f_Q_tktd(Q, A) * f_Dw_probit(Dw) - m_max - D) * A;
  }

  // Internal concentration of P
  dQ = v_max * f_Q_P_tktd(A, Q, P) * A - (m_max + D) * Q;

  // External concentration of P
  dP = D * (R_0 - P) + Q * m_max - (v_max * f_Q_P_tktd(A, Q, P) * A);

  // TKTD damage concentration
  dDw = kD * (Cw - Dw);

  // Additional outputs
  if(*ip > 0)
  {
    // external concentration
    yout[0] = Cw;
    // response functions
    if(*ip > 1) yout[1] = f_temp_tktd(T_act);
    if(*ip > 2) yout[2] = f_I_tktd(I);
    if(*ip > 3) yout[3] = f_Q_tktd(Q, A);
    if(*ip > 4) yout[4] = f_Q_P_tktd(A, Q, P);
    if(*ip > 5) yout[5] = dose_resp == 0 ? f_Dw_logit(Dw) : f_Dw_probit(Dw);
    // derivatives
    if(*ip > 6) yout[6] = dA;
    if(*ip > 7) yout[7] = dQ;
    if(*ip > 8) yout[8] = dP;
    if(*ip > 9) yout[9] = dDw;
  }
}
