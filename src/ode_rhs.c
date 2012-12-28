#include <gsl/gsl_errno.h>
#include "rate_coeffs.h"
#include "param.h"

/* Right-hand side of each ODE, i.e., the side with all the rates and
 * abundances multiplied together. See main.cpp for the isotope codes
 * since I've hard-coded them here. */

/* Function arguments are:
 * t = time (independent variable)
 * y[] = vector containing isotope abundances at time t
 * dydt[] = RHS of each ODE
 * params -> all parameters other than time (just temperature in
 *           this case) */

int
ode_rhs (double t, const double y[], double dydt[], void *params_in)
{

  // get the temperature
  struct param params = *(struct param *) params_in;
  const double T = params.T;

  dydt[0] =
    lambda_ijT_avg (6, 12, T, 'a') * y[6] * y[12] + lambda_ijT_avg (9, 12, T,
								    'a') *
    y[9] * y[12] + lambda_ijT (11, 12, T) * y[11] * y[12];
  dydt[1] =
    -lambda_ijT (1, 12, T) * y[1] * y[12] + lambda_ijT_avg (6, 12, T,
							    'a') * y[6] *
    y[12];
  dydt[2] = lambda_ijT (1, 12, T) * y[1] * y[12] - lambda_ij_beta (2) * y[2];
  dydt[3] = lambda_ij_beta (2) * y[2] - lambda_ijT (3, 12, T) * y[3] * y[12];
  dydt[4] =
    -lambda_ijT (4, 12, T) * y[4] * y[12] + lambda_ijT (3, 12,
							T) * y[3] * y[12] +
    lambda_ijT_avg (9, 12, T, 'a') * y[9] * y[12];
  dydt[5] = -lambda_ij_beta (5) * y[5] + lambda_ijT (4, 12, T) * y[4] * y[12];
  dydt[6] =
    lambda_ij_beta (5) * y[5] - lambda_ijT_avg (6, 12, T,
						'a') * y[6] * y[12] -
    lambda_ijT_avg (6, 12, T, 'g') * y[6] * y[12] + lambda_ijT (11, 12,
								T) * y[11] *
    y[12];
  dydt[7] =
    lambda_ijT_avg (6, 12, T, 'g') * y[6] * y[12] - lambda_ijT (7, 12,
								T) * y[7] *
    y[12];
  dydt[8] = lambda_ijT (7, 12, T) * y[7] * y[12] - lambda_ij_beta (8) * y[8];
  dydt[9] =
    lambda_ij_beta (8) * y[8] - lambda_ijT_avg (9, 12, T,
						'g') * y[9] * y[12] -
    lambda_ijT_avg (9, 12, T, 'a') * y[9] * y[12];
  dydt[10] =
    lambda_ijT_avg (9, 12, T,
		    'g') * y[9] * y[12] - lambda_ij_beta (10) * y[10];
  dydt[11] =
    lambda_ij_beta (10) * y[10] - lambda_ijT (11, 12, T) * y[11] * y[12];
  dydt[12] =
    -lambda_ijT (1, 12, T) * y[1] * y[12] - lambda_ijT (3, 12,
							T) * y[3] * y[12] -
    lambda_ijT (4, 12, T) * y[4] * y[12] - lambda_ijT_avg (6, 12, T,
							   'a') * y[6] *
    y[12] - lambda_ijT_avg (6, 12, T, 'g') * y[6] * y[12] - lambda_ijT (7, 12,
									T) *
    y[7] * y[12] - lambda_ijT_avg (9, 12, T,
				   'a') * y[9] * y[12] - lambda_ijT_avg (9,
									 12,
									 T,
									 'g')
    * y[9] * y[12] - lambda_ijT (11, 12, T) * y[11] * y[12];
  return GSL_SUCCESS;
}
