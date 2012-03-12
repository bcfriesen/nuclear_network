#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "rate_coeffs.hpp"
#include "jacobian.hpp"
#include "global.hpp"

/* Create sparse Jacobian matrix. For this CNO network the matrix is
 * 72% sparse.  Clearly a sparse solver would be in order here. I
 * tried SuperLU but it's impossible to use and can't do simple things
 * like adding two matrices together. Soooo, yeah, eff that. For small
 * jacobians (this only has 13^2 = 169 elements) the
 * performance/memory difference between sparse and full solvers is
 * immeasurably small. */

int jacobian (double t, const double y[], double *dfdy, double dfdt[],
              void *params) {

  // the only parameter is temperature
  double T = *(double *)params;
  /* a clever way to transform a regular 2-dimensional C array to a
   * GSL matrix */
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, nvar, nvar);
  gsl_matrix *jac = &dfdy_mat.matrix;

  // 1-D C array to GSL vector
  gsl_vector_const_view y_vec = gsl_vector_const_view_array(y, nvar);
  const gsl_vector *Y_i = &y_vec.vector;

  /* for now I'll ignore the time dependence in the beta-decay rates,
   * in which case the whole RHS has no explicit time dependence and
   * these derivatives are all zero. this shouldn't make any
   * difference since the rates are so insanely fast compared to the
   * 2-body rates. */
  for (unsigned int i = 0; i < nvar; ++i) { dfdt[i] = 0.0; }

  /* In my notes I started matrix elements at 1, not 0. So these
   * values are universally lower by 1. */

  /* the function gsl_matrix_set(mtx, i, j, val) just means "set the
   * matrix element mtx[i,j] = val". GSL matrices (and vectors) are
   * nice because they check for out-of-bounds calls, whereas regular
   * C arrays don't.
   *
   * similarly, the function gsl_vector_get(vec, i) means "give me the
   * value of vec[i]" */
  gsl_matrix_set(jac,  1,  1, -lambda_ij( 1, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  2,  1,  lambda_ij( 1, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 12,  1, -lambda_ij( 2, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  2,  2, -lambda_ij( 2));
  gsl_matrix_set(jac,  3,  2,  lambda_ij( 2));
  gsl_matrix_set(jac,  3,  3, -lambda_ij( 3, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  4,  3,  lambda_ij( 3, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 12,  3, -lambda_ij( 3, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  4,  4, -lambda_ij( 4, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  5,  4,  lambda_ij( 4, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 12,  4, -lambda_ij( 4, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  5,  5, -lambda_ij( 5));
  gsl_matrix_set(jac,  6,  5,  lambda_ij( 5));
  gsl_matrix_set(jac,  0,  6,  lambda_ij( 6, 12, T, 'a')
                             * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  1,  6,  lambda_ij( 6, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  6,  6, -lambda_ij( 6, 12, T, 'a')
                             * gsl_vector_get(Y_i, 12)
                             - lambda_ij( 6, 12, T, 'g')
			     * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  7,  6,  lambda_ij( 6, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 12,  6, -lambda_ij( 6, 12, T, 'a')
                             * gsl_vector_get(Y_i, 12)
                             - lambda_ij( 6, 12, T, 'g')
			     * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  7,  7, -lambda_ij( 7, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  8,  7,  lambda_ij( 7, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 12,  7, -lambda_ij( 7, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  8,  8, -lambda_ij( 8));
  gsl_matrix_set(jac,  9,  8,  lambda_ij( 8));
  gsl_matrix_set(jac,  0,  9,  lambda_ij( 9, 12, T, 'a')
                             * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  4,  9,  lambda_ij( 9, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  9,  9, -lambda_ij( 9, 12, T, 'g')
                             * gsl_vector_get(Y_i, 12)
                             - lambda_ij( 9, 12, T, 'a')
			     * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 10,  9,  lambda_ij( 9, 12, T, 'g')
                             * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 12,  9, -lambda_ij( 9, 12, T, 'a')
                             * gsl_vector_get(Y_i, 12)
                             - lambda_ij( 9, 12, T, 'g')
			     * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 10, 10, -lambda_ij(10));
  gsl_matrix_set(jac, 11, 10,  lambda_ij(10));
  gsl_matrix_set(jac,  0, 11,  lambda_ij(11, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  6, 11,  lambda_ij(11, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 11, 11, -lambda_ij(11, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac, 12, 11, -lambda_ij(11, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  0, 12,  lambda_ij( 9, 12, T,'a')
                             * gsl_vector_get(Y_i,  6)
                             + lambda_ij( 9, 12, T,'a')
                             * gsl_vector_get(Y_i,  9)
                             + lambda_ij(11, 12, T) * gsl_vector_get(Y_i, 11));
  gsl_matrix_set(jac,  1, 12, -lambda_ij( 1, 12, T) * gsl_vector_get(Y_i,  1)
                             + lambda_ij( 6, 12, T) * gsl_vector_get(Y_i,  6));
  gsl_matrix_set(jac,  2, 12,  lambda_ij( 1, 12, T) * gsl_vector_get(Y_i,  1));
  gsl_matrix_set(jac,  3, 12, -lambda_ij( 3, 12, T) * gsl_vector_get(Y_i,  3));
  gsl_matrix_set(jac,  4, 12, -lambda_ij( 4, 12, T) * gsl_vector_get(Y_i,  4)
                             + lambda_ij( 3, 12, T) * gsl_vector_get(Y_i,  3)
                             + lambda_ij( 9, 12, T) * gsl_vector_get(Y_i,  9));
  gsl_matrix_set(jac,  5, 12,  lambda_ij( 4, 12, T) * gsl_vector_get(Y_i,  4));
  gsl_matrix_set(jac,  6, 12, -lambda_ij( 6, 12, T, 'a')
                             * gsl_vector_get(Y_i,  6)
                             - lambda_ij( 6, 12, T, 'g')
			     * gsl_vector_get(Y_i,  6)
                             + lambda_ij(11, 12, T) * gsl_vector_get(Y_i, 12));
  gsl_matrix_set(jac,  7, 12,  lambda_ij( 6, 12, T) * gsl_vector_get(Y_i,  6)
                             - lambda_ij( 7, 12, T) * gsl_vector_get(Y_i,  7));
  gsl_matrix_set(jac,  8, 12,  lambda_ij( 7, 12, T) * gsl_vector_get(Y_i,  7));
  gsl_matrix_set(jac,  9, 12, -lambda_ij( 9, 12, T, 'g')
                             * gsl_vector_get(Y_i,  9)
                             - lambda_ij( 9, 12, T, 'a')
			     * gsl_vector_get(Y_i,  9));
  gsl_matrix_set(jac, 10, 12,  lambda_ij( 9, 12, T, 'g')
                             * gsl_vector_get(Y_i,  9));
  gsl_matrix_set(jac, 11, 12, -lambda_ij(11, 12, T)
                             * gsl_vector_get(Y_i, 11));
  gsl_matrix_set(jac, 12, 12, -lambda_ij( 1, 12, T) * gsl_vector_get(Y_i,  1)
                             - lambda_ij( 3, 12, T) * gsl_vector_get(Y_i,  3)
                             - lambda_ij( 4, 12, T) * gsl_vector_get(Y_i,  4)
		             - lambda_ij( 6, 12, T, 'a')
			     * gsl_vector_get(Y_i,  6)
	                     - lambda_ij( 6, 12, T, 'g')
			     * gsl_vector_get(Y_i,  6)
		             - lambda_ij( 7, 12, T) * gsl_vector_get(Y_i,  7)
	                     - lambda_ij( 9, 12, T, 'g')
			     * gsl_vector_get(Y_i,  9)
		             - lambda_ij( 9, 12, T, 'a')
			     * gsl_vector_get(Y_i,  9)
	                     - lambda_ij( 11, 12, T) * gsl_vector_get(Y_i, 11));
  return GSL_SUCCESS;
}
