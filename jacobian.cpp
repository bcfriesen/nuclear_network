#include <gsl/gsl_matrix.h>
#include "rate_coeffs.hpp"
#include "jacobian.hpp"

#define N_NONZERO_EL 47
#define N_ISO 13
/* Create sparse Jacobian matrix. For this CNO network the matrix is 72% sparse.
 * Clearly a sparse solver would be in order here. I tried SuperLU but it's
 * impossible to use and can't do simple things like adding two matrices
 * together. Eff that. */

int create_jacobian_matrix (gsl_matrix *jac, gsl_vector *Y_i, double T) {
  /* In my notes I started matrix elements at 1, not 0. So these values are
   * universally lower by 1. */

  // J_{1, 1}
  gsl_matrix_set(jac, 1, 1, -lambda_ij(1, 12, T) * gsl_vector_get(Y_i, 12));
  // J_{2, 1}
  gsl_matrix_set(jac, 2, 1, lambda_ij(1, 12, T) * gsl_vector_get(Y_i, 12));
  // J_{12, 1}
  gsl_matrix_set(jac, 12, 1, -lambda_ij(1, 12, T) * gsl_vector_get(Y_i, 12));
  // J_{2, 2}
  gsl_matrix_set(jac, 2, 2, -lambda_ij(2));
  // J_{3, 2}
  gsl_matrix_set(jac, 3, 2, lambda_ij(2));
  // J_{3, 3}
  gsl_matrix_set(jac, 3, 3, -lambda_ij(3, 12, T) * gsl_vector_get(Y_i, 12));
  // J_{4, 3}
  gsl_matrix_set(jac, 4, 3, lambda_ij(3, 12, T) * gsl_vector_get(Y_i, 12));
  // J_{12, 3}
  gsl_matrix_set(jac, 12, 3, -lambda_ij(3,12,T) * gsl_vector_get(Y_i, 12));
  // J_{4, 4}
  gsl_matrix_set(jac, 4, 4, -lambda_ij(4,12,T) * gsl_vector_get(Y_i, 12));
  // J_{5, 4}
  gsl_matrix_set(jac, 5, 4, lambda_ij(4,12,T) * gsl_vector_get(Y_i, 12));
  // J_{12, 4}
  gsl_matrix_set(jac, 12, 4, -lambda_ij(4,12,T) * gsl_vector_get(Y_i, 12));
  // J_{5, 5}
  gsl_matrix_set(jac, 5, 5, -lambda_ij(5));
  // J_{6, 5}
  gsl_matrix_set(jac, 6, 5, lambda_ij(5));
  // J_{0, 6}
  gsl_matrix_set(jac, 0, 6, lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 12));
  // J_{1, 6}
  gsl_matrix_set(jac, 1, 6, lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 12));
  // J_{6, 6}
  gsl_matrix_set(jac, 6, 6, -lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 12)
                 - lambda_ij(6,12,T,'g') * gsl_vector_get(Y_i, 12));
  // J_{7, 6}
  gsl_matrix_set(jac, 7, 6, lambda_ij(6,12,T,'g') * gsl_vector_get(Y_i, 12));
  // J_{12, 6}
  gsl_matrix_set(jac, 12, 6, -lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 12)
                 -  lambda_ij(6,12,T,'g') * gsl_vector_get(Y_i, 12));
  // J_{7, 7}
  gsl_matrix_set(jac, 7, 7, -lambda_ij(7,12,T) * gsl_vector_get(Y_i, 12));
  // J_{8, 7}
  gsl_matrix_set(jac, 8, 7, lambda_ij(7,12,T) * gsl_vector_get(Y_i, 12));
  // J_{12, 7}
  gsl_matrix_set(jac, 12, 7, -lambda_ij(7,12,T) * gsl_vector_get(Y_i, 12));
  // J_{8, 8}
  gsl_matrix_set(jac, 8, 8, -lambda_ij(8));
  // J_{9, 8}
  gsl_matrix_set(jac, 9, 8, lambda_ij(8));
  // J_{0, 9}
  gsl_matrix_set(jac, 0, 9, lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 12));
  // J_{4, 9}
  gsl_matrix_set(jac, 4, 9, lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 12));
  // J_{9, 9}
  gsl_matrix_set(jac, 9, 9, -lambda_ij(9,12,T,'g') * gsl_vector_get(Y_i, 12)
                 -  lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 12));
  // J_{10, 9}
  gsl_matrix_set(jac, 10, 9, lambda_ij(9,12,T,'g') * gsl_vector_get(Y_i, 12));
  // J_{12, 9}
  gsl_matrix_set(jac, 12, 9, -lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 12)
                 -  lambda_ij(9,12,T,'g') * gsl_vector_get(Y_i, 12));
  // J_{10, 10}
  gsl_matrix_set(jac, 10, 10, -lambda_ij(10));
  // J_{11, 10}
  gsl_matrix_set(jac, 11, 10, lambda_ij(10));
  // J_{0, 11}
  gsl_matrix_set(jac, 0, 11, lambda_ij(11,12,T) * gsl_vector_get(Y_i, 12));
  // J_{6, 11}
  gsl_matrix_set(jac, 6, 11, lambda_ij(11,12,T) * gsl_vector_get(Y_i, 12));
  // J_{11, 11}
  gsl_matrix_set(jac, 11, 11, -lambda_ij(11,12,T) * gsl_vector_get(Y_i, 12));
  // J_{12, 11}
  gsl_matrix_set(jac, 12, 11, -lambda_ij(11,12,T) * gsl_vector_get(Y_i, 12));
  // J_{0, 12}
  gsl_matrix_set(jac, 0, 12,  lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 6)
                 + lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 9)
                 + lambda_ij(11,12,T) * gsl_vector_get(Y_i, 11));
  // J_{1, 12}
  gsl_matrix_set(jac, 1, 12, -lambda_ij(1,12,T) * gsl_vector_get(Y_i, 1)
                 + lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 6));
  // J_{2, 12}
  gsl_matrix_set(jac, 2, 12, lambda_ij(1,12,T) * gsl_vector_get(Y_i, 1));
  // J_{3, 12}
  gsl_matrix_set(jac, 3, 12, -lambda_ij(3,12,T) * gsl_vector_get(Y_i, 3));
  // J_{4, 12}
  gsl_matrix_set(jac, 4, 12, -lambda_ij(4,12,T) * gsl_vector_get(Y_i, 4)
                 + lambda_ij(3,12,T) * gsl_vector_get(Y_i, 3)
                 + lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 9));
  // J_{5, 12}
  gsl_matrix_set(jac, 5, 12, lambda_ij(4,12,T) * gsl_vector_get(Y_i, 4));
  // J_{6, 12}
  gsl_matrix_set(jac, 6, 12, -lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 6)
                 - lambda_ij(6,12,T,'g') * gsl_vector_get(Y_i, 6)
                 + lambda_ij(11,12,T) * gsl_vector_get(Y_i, 12));
  // J_{7, 12}
  gsl_matrix_set(jac, 7, 12, lambda_ij(6,12,T,'g') * gsl_vector_get(Y_i, 6)
                 - lambda_ij(7,12,T) * gsl_vector_get(Y_i, 7));
  // J_{8, 12}
  gsl_matrix_set(jac, 8, 12, lambda_ij(7,12,T) * gsl_vector_get(Y_i, 7));
  // J_{9, 12}
  gsl_matrix_set(jac, 9, 12, -lambda_ij(9,12,T,'g') * gsl_vector_get(Y_i, 9)
                 - lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 9));
  // J_{10, 12}
  gsl_matrix_set(jac, 10, 12, lambda_ij(9,12,T,'g') * gsl_vector_get(Y_i, 9));
  // J_{11, 12}
  gsl_matrix_set(jac, 11, 12, -lambda_ij(11,12,T) * gsl_vector_get(Y_i, 12));
  // J_{12, 12}
  gsl_matrix_set(jac, 12, 12, -lambda_ij(1,12,T) * gsl_vector_get(Y_i, 1)
                 - lambda_ij(3,12,T) * gsl_vector_get(Y_i, 3)
                 - lambda_ij(4,12,T) * gsl_vector_get(Y_i, 4)
		 - lambda_ij(6,12,T,'a') * gsl_vector_get(Y_i, 6)
	         - lambda_ij(6,12,T,'g') * gsl_vector_get(Y_i, 6)
		 - lambda_ij(7,12,T) * gsl_vector_get(Y_i, 7)
	         - lambda_ij(9,12,T,'g') * gsl_vector_get(Y_i, 9)
		 - lambda_ij(9,12,T,'a') * gsl_vector_get(Y_i, 9)
	         - lambda_ij(11,12,T) * gsl_vector_get(Y_i, 11));
  return 0;
}
