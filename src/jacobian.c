#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "rate_coeffs.h"
#include "global.h"

/* Create sparse Jacobian matrix. For this CNO network the matrix is 72%
 * sparse. A sparse solver would be in order here; I tried SuperLU but it's
 * impossible to use and can't do simple things like adding two matrices
 * together. For small jacobians (this only has 13^2 = 169 elements) the
 * performance/memory difference between sparse and full solvers is small. */

/* the arguments of this function are:
 * t -> time (independent variable)
 * y[] -> vector containing relative abundances (by number) of each isotope at
 *        time t
 * dfdy -> Jacobian matrix dfdt -> partial derivative of RHS of each ODE w.r.t.
 *         time
 * params -> any arguments that the Jacobian matrix elements may need
 * besides the independent variable (time). In this case the only auxiliary
 * parameter is temperature */
int jacobian (double t, const double y[], double *dfdy, double dfdt[],
              void *params) {
  // loops
  unsigned int i;
  // the only parameter is temperature
  double T = *(double *)params;

  /* for now I'll ignore the time dependence in the beta-decay rates,
   * in which case the whole RHS has no explicit time dependence and
   * these derivatives are all zero. this shouldn't make any
   * difference since the rates are so insanely fast compared to the
   * 2-body rates. */
  for (i = 0; i < nvar; ++i) { dfdt[i] = 0.0; }

  /* GSL expects the Jacobian matrix to be stored in row-major order in a 1-D
   * vector, so J[i][j] = dfdy[i*DIM + j]. Hence the weird notation here. */
  dfdy[ 1*nvar +  1] = -lambda_ijT( 1, 12, T) * y[12];
  dfdy[ 2*nvar +  1] =  lambda_ijT( 1, 12, T) * y[12];
  dfdy[12*nvar +  1] = -lambda_ijT( 2, 12, T) * y[12];
  dfdy[ 2*nvar +  2] = -lambda_ij_beta( 2);
  dfdy[ 3*nvar +  2] =  lambda_ij_beta( 2);
  dfdy[ 3*nvar +  3] = -lambda_ijT( 3, 12, T) * y[12];
  dfdy[ 4*nvar +  3] =  lambda_ijT( 3, 12, T) * y[12];
  dfdy[12*nvar +  3] = -lambda_ijT( 3, 12, T) * y[12];
  dfdy[ 4*nvar +  4] = -lambda_ijT( 4, 12, T) * y[12];
  dfdy[ 5*nvar +  4] =  lambda_ijT( 4, 12, T) * y[12];
  dfdy[12*nvar +  4] = -lambda_ijT( 4, 12, T) * y[12];
  dfdy[ 5*nvar +  5] = -lambda_ij_beta( 5);
  dfdy[ 6*nvar +  5] =  lambda_ij_beta( 5);
  dfdy[ 0*nvar +  6] =  lambda_ijT_avg( 6, 12, T, 'a') * y[12];
  dfdy[ 1*nvar +  6] =  lambda_ijT( 6, 12, T) * y[12];
  dfdy[ 6*nvar +  6] = -lambda_ijT_avg( 6, 12, T, 'a') * y[12] - lambda_ijT_avg( 6, 12, T, 'g') * y[12];
  dfdy[ 7*nvar +  6] =  lambda_ijT( 6, 12, T) * y[12];
  dfdy[12*nvar +  6] = -lambda_ijT_avg( 6, 12, T, 'a') * y[12] - lambda_ijT_avg( 6, 12, T, 'g') * y[12];
  dfdy[ 7*nvar +  7] = -lambda_ijT( 7, 12, T) * y[12];
  dfdy[ 8*nvar +  7] =  lambda_ijT( 7, 12, T) * y[12];
  dfdy[12*nvar +  7] = -lambda_ijT( 7, 12, T) * y[12];
  dfdy[ 8*nvar +  8] = -lambda_ij_beta( 8);
  dfdy[ 9*nvar +  8] =  lambda_ij_beta( 8);
  dfdy[ 0*nvar +  9] =  lambda_ijT_avg( 9, 12, T, 'a') * y[12];
  dfdy[ 4*nvar +  9] =  lambda_ijT( 9, 12, T) * y[12];
  dfdy[ 9*nvar +  9] = -lambda_ijT_avg( 9, 12, T, 'g') * y[12] - lambda_ijT_avg( 9, 12, T, 'a') * y[12];
  dfdy[10*nvar +  9] =  lambda_ijT_avg( 9, 12, T, 'g') * y[12];
  dfdy[12*nvar +  9] = -lambda_ijT_avg( 9, 12, T, 'a') * y[12] - lambda_ijT_avg( 9, 12, T, 'g') * y[12];
  dfdy[10*nvar + 10] = -lambda_ij_beta(10);
  dfdy[11*nvar + 10] =  lambda_ij_beta(10);
  dfdy[ 0*nvar + 11] =  lambda_ijT(11, 12, T) * y[12];
  dfdy[ 6*nvar + 11] =  lambda_ijT(11, 12, T) * y[12];
  dfdy[11*nvar + 11] = -lambda_ijT(11, 12, T) * y[12];
  dfdy[12*nvar + 11] = -lambda_ijT(11, 12, T) * y[12];
  dfdy[ 0*nvar + 12] =  lambda_ijT_avg( 9, 12, T,'a') * y[6] + lambda_ijT_avg( 9, 12, T,'a') * y[9] + lambda_ijT(11, 12, T) * y[11];
  dfdy[ 1*nvar + 12] = -lambda_ijT( 1, 12, T) * y[1] + lambda_ijT( 6, 12, T) * y[6];
  dfdy[ 2*nvar + 12] =  lambda_ijT( 1, 12, T) * y[1];
  dfdy[ 3*nvar + 12] = -lambda_ijT( 3, 12, T) * y[3];
  dfdy[ 4*nvar + 12] = -lambda_ijT( 4, 12, T) * y[4] + lambda_ijT( 3, 12, T) * y[3] + lambda_ijT( 9, 12, T) * y[9];
  dfdy[ 5*nvar + 12] =  lambda_ijT( 4, 12, T) * y[4];
  dfdy[ 6*nvar + 12] = -lambda_ijT_avg( 6, 12, T, 'a') * y[6] - lambda_ijT_avg( 6, 12, T, 'g') * y[6] + lambda_ijT(11, 12, T) * y[12];
  dfdy[ 7*nvar + 12] =  lambda_ijT( 6, 12, T) * y[6] - lambda_ijT( 7, 12, T) * y[7];
  dfdy[ 8*nvar + 12] =  lambda_ijT( 7, 12, T) * y[7];
  dfdy[ 9*nvar + 12] = -lambda_ijT_avg( 9, 12, T, 'g') * y[9] - lambda_ijT_avg( 9, 12, T, 'a') * y[9];
  dfdy[10*nvar + 12] =  lambda_ijT_avg( 9, 12, T, 'g') * y[9];
  dfdy[11*nvar + 12] = -lambda_ijT(11, 12, T) * y[11];
  dfdy[12*nvar + 12] = -lambda_ijT( 1, 12, T) * y[1] - lambda_ijT( 3, 12, T) * y[3] - lambda_ijT( 4, 12, T) * y[4] - lambda_ijT_avg( 6, 12, T, 'a') * y[6] - lambda_ijT_avg( 6, 12, T, 'g') * y[6] - lambda_ijT( 7, 12, T) * y[7] - lambda_ijT_avg( 9, 12, T, 'g') * y[9] - lambda_ijT_avg( 9, 12, T, 'a') * y[9] - lambda_ijT( 11, 12, T) * y[11];
  return GSL_SUCCESS;
}
