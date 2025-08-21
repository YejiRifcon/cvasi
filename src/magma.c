#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>
/*******************************************************************
 *
 * The Magma model (Witt et al., submitted) interprets the Tier 2C
 * version of the Lemna model by Klein et al. (2021).
 *
 * The model provides additional output on intermediary variables on
 * request; please refer to deSolve's manual on the 'nout' parameter.
 *
 *******************************************************************/

/**
 * Allocate memory for global parameter array
 */
static double parms[13] = {0};

/**
 * Allocate memory for forcing function data
 *
 * Array's' values get updated by ODE solver in every time step.
 */
static double forc[1] = {0};
// Internal variable for growth model selection
static short growth_model = 0; // magic value: 0, invalid dummy value

/*
 * Define aliases
 */
// state variable aliases
#define BM    y[0]
#define M_int y[1]

// derivative aliases
#define dBM    ydot[0]
#define dM_int ydot[1]

// parameter aliases
// growth model parameters
#define mu_control   parms[0]
#define growthno     parms[1]
#define D_L          parms[2]
// toxicodynamic parameters
#define E_max         parms[3]
#define EC50_int      parms[4]
#define b             parms[5]
// toxicokinetic parameters
#define P             parms[6]
#define r_A_DW        parms[7]
#define r_FW_DW       parms[8]
#define r_FW_V        parms[9]
#define r_DW_TSL      parms[10]
#define K_pw          parms[11]
#define k_met         parms[12]
// forcings by environmental variables
#define C_ext forc[0]
// constants
#define EXP_GROWTH 1 // magic value: constant set in R solver code
#define LOG_GROWTH 2 // magic value: constant set in R solver code

/**
 * Parameter initializer
 */
void magma_init(void (* odeparms)(int *, double *))
{
  int N=13;
  odeparms(&N, parms);

  growth_model = (short)growthno;
}

/**
 * Forcings initializer
 */
void magma_forc(void (* odeforcs)(int *, double *))
{
  int N=1;
  odeforcs(&N, forc);
}

/* Concentration response of photosynthesis [Toxicodynamics] (Box 9)
 * @param C_int_unb internal unbound toxicant concentration (mass per volume, e.g. ug L-1)
 * @param E_max maximum inhibition (-)
 * @param EC50_int int. conc. resulting in 50% effect (mass per volume, e.g. ug L-1)
 * @param b slope parameter (-)
 * @return value from the interval [0,1]
 */
double fCint_photo_(double C_int_unb) {
  double pow_C_int_b = pow(C_int_unb, b);
  return(1 - E_max * pow_C_int_b / (pow(EC50_int, b) + pow_C_int_b));
}

/**
 * Derivatives
 */
void magma_func(int *neq, double *t, double *y, double *ydot, double *yout, int*ip)
{
  //
  // Toxicokinetics
  //

  // Internal toxicant concentration (ug L-1)
  double C_int, C_int_unb;
  if(BM <= 0) { // avoid division by zero
    C_int = 0;
    C_int_unb = 0;
  } else {
    C_int = M_int * r_FW_V / (BM * r_FW_DW);
    C_int_unb = C_int / K_pw; // unbound internal concentration
  }

  // TK model ODE
  dM_int = P * BM * r_A_DW * (C_ext - C_int_unb) - M_int / K_pw * k_met;

  //
  // Effects on photosynthesis
  //

  // Photosynthesis dependency function
  double f_photo = fCint_photo_(C_int_unb);

  //
  // Population growth
  //

  // Growth model ODE
  if(growth_model == EXP_GROWTH) {
    dBM = mu_control * f_photo * BM; // exponential growth
  } else if(growth_model == LOG_GROWTH) {
    dBM = mu_control * f_photo * BM * (1 - BM / D_L); // logistic growth
  } else {
    Rf_error("unknown growth function selected");
  }

  //
  // Additional output variables
  //
  if(*ip > 0)
  {
    // total shoot length
    yout[0] = BM / r_DW_TSL;
    // internal concentration
    if(*ip > 1) yout[1] = C_int;
    // response function
    if(*ip > 2) yout[2] = f_photo;
    // toxicant concentration
    if(*ip > 3) yout[3] = C_int_unb;
    if(*ip > 4) yout[4] = C_ext;
    // derivatives
    if(*ip > 5) yout[5] = dBM;
    if(*ip > 6) yout[6] = dM_int;
  }
}
