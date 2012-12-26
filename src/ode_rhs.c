#include <gsl/gsl_matrix.h>
#include "rate_coeffs.h"
#include "global.h"

/* Right-hand side of each ODE, i.e., the side with all the rates and
 * abundances multiplied together. See main.cpp for the isotope codes
 * since I've hard-coded them here. */

/* Function arguments are:
 * t = time (independent variable)
 * y[] = vector containing isotope abundances at time t
 * dydt[] = RHS of each ODE
 * params -> all parameters other than time (just temperature in
 *           this case) */

int ode_rhs(double t, const double y[], double dydt[], void *params) {
  // convert C arrays to GSL vectors
  gsl_vector_view dydt_vec = gsl_vector_view_array(dydt, nvar);
  gsl_vector *ode_rhs = &dydt_vec.vector;

  // same here
  gsl_vector_const_view y_vec = gsl_vector_const_view_array(y, nvar);
  const gsl_vector *Y_i = &y_vec.vector;

  // get the temperature
  double T = *(double *)params;

  // the RHS of each ODE
  gsl_vector_set(ode_rhs, 0, lambda_ijT_avg(6, 12, T, 'a') * gsl_vector_get(Y_i, 6)
                 * gsl_vector_get(Y_i, 12)
		 + lambda_ijT_avg(9, 12, T, 'a') * gsl_vector_get(Y_i, 9)
		 * gsl_vector_get(Y_i, 12)
		 + lambda_ijT(11, 12, T) * gsl_vector_get(Y_i, 11)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 1, -lambda_ijT(1, 12, T) * gsl_vector_get(Y_i, 1)
                 * gsl_vector_get(Y_i, 12)
		 + lambda_ijT_avg(6, 12, T, 'a') * gsl_vector_get(Y_i, 6)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 2, lambda_ijT(1, 12, T) * gsl_vector_get(Y_i, 1)
                 * gsl_vector_get(Y_i, 12)
		 - lambda_ij_beta(2) * gsl_vector_get(Y_i, 2));

  gsl_vector_set(ode_rhs, 3, lambda_ij_beta(2) * gsl_vector_get(Y_i, 2)
                 - lambda_ijT(3, 12, T) * gsl_vector_get(Y_i, 3)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 4, -lambda_ijT(4, 12, T) * gsl_vector_get(Y_i, 4)
                 * gsl_vector_get(Y_i, 12)
		 + lambda_ijT(3, 12, T) * gsl_vector_get(Y_i, 3)
		 * gsl_vector_get(Y_i, 12)
		 + lambda_ijT_avg(9, 12, T, 'a') * gsl_vector_get(Y_i, 9)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 5, -lambda_ij_beta(5) * gsl_vector_get(Y_i, 5)
                 + lambda_ijT(4, 12, T) * gsl_vector_get(Y_i, 4)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 6, lambda_ij_beta(5) * gsl_vector_get(Y_i, 5)
                 - lambda_ijT_avg(6, 12, T, 'a') * gsl_vector_get(Y_i, 6)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT_avg(6, 12, T, 'g') * gsl_vector_get(Y_i, 6)
		 * gsl_vector_get(Y_i, 12)
		 + lambda_ijT(11, 12, T) * gsl_vector_get(Y_i, 11)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 7, lambda_ijT_avg(6, 12, T, 'g') * gsl_vector_get(Y_i, 6)
                 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT(7, 12, T) * gsl_vector_get(Y_i, 7)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 8, lambda_ijT(7, 12, T) * gsl_vector_get(Y_i, 7)
                 * gsl_vector_get(Y_i, 12)
		 - lambda_ij_beta(8) * gsl_vector_get(Y_i, 8));

  gsl_vector_set(ode_rhs, 9, lambda_ij_beta(8) * gsl_vector_get(Y_i, 8)
                 - lambda_ijT_avg(9, 12, T, 'g') * gsl_vector_get(Y_i, 9)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT_avg(9, 12, T, 'a') * gsl_vector_get(Y_i, 9)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 10, lambda_ijT_avg(9, 12, T, 'g') * gsl_vector_get(Y_i, 9)
                 * gsl_vector_get(Y_i, 12)
		 - lambda_ij_beta(10) * gsl_vector_get(Y_i, 10));

  gsl_vector_set(ode_rhs, 11, lambda_ij_beta(10) * gsl_vector_get(Y_i, 10)
                 - lambda_ijT(11, 12, T) * gsl_vector_get(Y_i, 11)
		 * gsl_vector_get(Y_i, 12));

  gsl_vector_set(ode_rhs, 12, -lambda_ijT(1, 12, T) * gsl_vector_get(Y_i, 1)
                 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT(3, 12, T) * gsl_vector_get(Y_i, 3)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT(4, 12, T) * gsl_vector_get(Y_i, 4)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT_avg(6, 12, T, 'a') * gsl_vector_get(Y_i, 6)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT_avg(6, 12, T, 'g') * gsl_vector_get(Y_i, 6)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT(7, 12, T) * gsl_vector_get(Y_i, 7)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT_avg(9, 12, T, 'a') * gsl_vector_get(Y_i, 9)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT_avg(9, 12, T, 'g') * gsl_vector_get(Y_i, 9)
		 * gsl_vector_get(Y_i, 12)
		 - lambda_ijT(11, 12, T) * gsl_vector_get(Y_i, 11)
		 * gsl_vector_get(Y_i, 12));
  return GSL_SUCCESS;
}
