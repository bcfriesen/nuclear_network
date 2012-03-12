#include "rate_coeffs.hpp"
#include "jacobian.hpp"

/* Create sparse Jacobian matrix. For this CNO network the matrix is 72% sparse.
 * Clearly a sparse solver would be in order here. I tried SuperLU but it's
 * impossible to use and can't do simple things like adding two matrices
 * together. Eff that. */

int jacobian (double t, const double y[], double *dfdy, double dfdt[],
              void *params) {
  double temp = *(double *)params;
  /* In my notes I started matrix elements at 1, not 0. So these values are
   * universally lower by 1. */

  dfdy[ 1][ 1] = -lambda_ij( 1, 12, temp     ) * y[12];
  dfdy[ 2][ 1] =  lambda_ij( 1, 12, temp     ) * y[12];
  dfdy[12][ 1] = -lambda_ij( 2, 12, temp     ) * y[12];
  dfdy[ 2][ 2] = -lambda_ij( 2               ) * y[ 2];
  dfdy[ 3][ 2] =  lambda_ij( 2               ) * y[ 2];
  dfdy[ 3][ 3] = -lambda_ij( 3, 12, temp     ) * y[12];
  dfdy[ 4][ 3] =  lambda_ij( 3, 12, temp     ) * y[12];
  dfdy[12][ 3] = -lambda_ij( 3, 12, temp     ) * y[12];
  dfdy[ 4][ 4] = -lambda_ij( 4, 12, temp     ) * y[12];
  dfdy[ 5][ 4] =  lambda_ij( 4, 12, temp     ) * y[12];
  dfdy[12][ 4] = -lambda_ij( 4, 12, temp     ) * y[12];
  dfdy[ 5][ 5] = -lambda_ij( 5               ) * y[ 5];
  dfdy[ 6][ 5] =  lambda_ij( 5               ) * y[ 5];
  dfdy[ 0][ 6] =  lambda_ij( 6, 12, temp, 'a') * y[12];
  dfdy[ 1][ 6] =  lambda_ij( 6, 12, temp     ) * y[12];
  dfdy[ 6][ 6] =  lambda_ij( 6, 12, temp, 'a') * y[12]
                - lambda_ij( 6, 12, temp, 'g') * y[12];
  dfdy[ 7][ 6] =  lambda_ij( 6, 12, temp     ) * y[12];
  dfdy[12][ 6] = -lambda_ij( 6, 12, temp, 'a') * y[12]
                - lambda_ij( 6, 12, temp, 'g') * y[12];
  dfdy[ 7][ 7] = -lambda_ij( 7, 12, temp     ) * y[12];
  dfdy[ 8][ 7] =  lambda_ij( 7, 12, temp     ) * y[12];
  dfdy[12][ 7] = -lambda_ij( 7, 12, temp     ) * y[12];
  dfdy[ 8][ 8] = -lambda_ij( 8               ) * y[ 8];
  dfdy[ 9][ 8] =  lambda_ij( 8               ) * y[ 8];
  dfdy[ 0][ 9] =  lambda_ij( 9, 12, temp, 'a') * y[12];
  dfdy[ 4][ 9] =  lambda_ij( 9, 12, temp     ) * y[12];
  dfdy[ 9][ 9] = -lambda_ij( 9, 12, temp, 'g') * y[12]
                - lambda_ij( 9, 12, temp, 'a') * y[12];
  dfdy[10][ 9] =  lambda_ij( 9, 12, temp, 'g') * y[12];
  dfdy[12][ 9] = -lambda_ij( 9, 12, temp, 'a') * y[12]
                - lambda_ij( 9, 12, temp, 'g') * y[12];
  dfdy[10][10] = -lambda_ij(10               ) * y[10];
  dfdy[11][10] =  lambda_ij(10               ) * y[10];
  dfdy[ 0][11] =  lambda_ij(11, 12, temp     ) * y[12];
  dfdy[ 6][11] =  lambda_ij(11, 12, temp     ) * y[12];
  dfdy[11][11] = -lambda_ij(11, 12, temp     ) * y[12];
  dfdy[12][11] = -lambda_ij(11, 12, temp     ) * y[12];
  dfdy[ 0][12] =  lambda_ij( 9, 12, temp, 'a') * y[ 6]
                + lambda_ij( 9, 12, temp, 'a') * y[ 9]
		+ lambda_ij(11, 12, temp     ) * y[11];
  dfdy[ 1][12] = -lambda_ij( 1, 12, temp     ) * y[ 1]
                + lambda_ij( 6, 12, temp     ) * y[ 6];
  dfdy[ 2][12] =  lambda_ij( 1, 12, temp     ) * y[ 1];
  dfdy[ 3][12] = -lambda_ij( 3, 12, temp     ) * y[ 3];
  dfdy[ 4][12] = -lambda_ij( 4, 12, temp     ) * y[ 4]
                + lambda_ij( 3, 12, temp     ) * y[ 3]
		+ lambda_ij( 9, 12, temp     ) * y[ 9];
  dfdy[ 5][12] =  lambda_ij( 4, 12, temp     ) * y[ 4];
  dfdy[ 6][12] = -lambda_ij( 6, 12, temp, 'a') * y[ 6]
                - lambda_ij( 6, 12, temp, 'g') * y[ 6]
		+ lambda_ij(11, 12, temp     ) * y[12];
  dfdy[ 7][12] =  lambda_ij( 6, 12, temp     ) * y[ 6]
                - lambda_ij( 7, 12, temp     ) * y[ 7];
  dfdy[ 8][12] =  lambda_ij( 7, 12, temp     ) * y[ 7];
  dfdy[ 9][12] = -lambda_ij( 9, 12, temp, 'g') * y[ 9]
                - lambda_ij( 9, 12, temp, 'a') * y[ 9];
  dfdy[10][12] =  lambda_ij( 9, 12, temp, 'g') * y[ 9];
  dfdy[11][12] = -lambda_ij(11, 12, temp     ) * y[11];
  dfdy[12][12] = -lambda_ij( 1, 12, temp     ) * y[ 1]
                - lambda_ij( 3, 12, temp     ) * y[ 3]
                - lambda_ij( 4, 12, temp     ) * y[ 4]
                - lambda_ij( 6, 12, temp, 'a') * y[ 6]
                - lambda_ij( 6, 12, temp, 'g') * y[ 6]
                - lambda_ij( 7, 12, temp     ) * y[ 7]
                - lambda_ij( 9, 12, temp, 'a') * y[ 9]
                - lambda_ij( 9, 12, temp, 'g') * y[ 9]
                - lambda_ij(11, 12, temp     ) * y[11]
  return 0;
}
