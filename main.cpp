#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include "jacobian.hpp"
#include "ode_rhs.hpp"
#include "rate_coeffs.hpp"

#define N_NONZERO_EL 47
#define N_ISO 13

int main() {
  const double T = 15.0e+06;

  unsigned int i;
  double temp;

  int s, error;

  // time step (sec)
  double h = 1.0e-7;
  // start time (sec)
  double t_now = 0.0;
  // stop time (sec)
  const double t_stop = 0.01;

  // declare a bunch of variables
  gsl_vector *A_i = gsl_vector_alloc(N_ISO);
  gsl_matrix *jac = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_matrix *lhs = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_vector *Y_i = gsl_vector_calloc(N_ISO);
  gsl_vector *ode_rhs = gsl_vector_calloc(N_ISO);
  gsl_vector *delta = gsl_vector_calloc(N_ISO);
  gsl_permutation *p = gsl_permutation_alloc(N_ISO);

  // set the atomic numbers (since Y_i = X_i / A_i, where Y_i is molar fraction,
  // X_i is mass fraction, and A_i is atomic #)
  gsl_vector_set(A_i, 0, 4.0);
  gsl_vector_set(A_i, 1,12.0);
  gsl_vector_set(A_i, 2,13.0);
  gsl_vector_set(A_i, 3,13.0);
  gsl_vector_set(A_i, 4,14.0);
  gsl_vector_set(A_i, 5,15.0);
  gsl_vector_set(A_i, 6,15.0);
  gsl_vector_set(A_i, 7,16.0);
  gsl_vector_set(A_i, 8,17.0);
  gsl_vector_set(A_i, 9,17.0);
  gsl_vector_set(A_i,10,18.0);
  gsl_vector_set(A_i,11,18.0);
  gsl_vector_set(A_i,12, 1.0);

  // mostly H1, with a tiny bit of C12, otherwise the rxn won't go!
  gsl_vector_set(Y_i, N_ISO-1, 0.99);
  gsl_vector_set(Y_i, 1, 8.3e-4);

  // print column headers
  printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s "
         "%12s\n",
         "t_now", "He4", "C12", "N13", "C13", "N14", "O15", "N15", "O16", "F17", 
	 "O17", "F18", "O18", "H1");
  do {
    error = create_jacobian_matrix(jac, Y_i, T);

    // create matrix to be inverted
    gsl_matrix_set_identity(lhs);
    for (unsigned int i = 0; i < N_ISO; i++) {
      gsl_matrix_set(lhs, i, i, gsl_matrix_get(lhs, i, i) / h);
    }
    error = gsl_matrix_sub(lhs, jac);

    // create RHS of matrix equation
    error = build_ode_rhs(ode_rhs, Y_i, T);

    // solve matrix equation
    gsl_linalg_LU_decomp(lhs, p, &s);
    gsl_linalg_LU_solve(lhs, p, ode_rhs, delta);

    // update abundances
    gsl_vector_add(Y_i, delta);

    // print abundances
    printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e "
           "%12.6e %12.6e %12.6e %12.6e %12.6e\n",
           t_now,
           gsl_vector_get(Y_i, 0)*4.0,
           gsl_vector_get(Y_i, 1)*12.0,
           gsl_vector_get(Y_i, 2)*13.0,
           gsl_vector_get(Y_i, 3)*13.0,
           gsl_vector_get(Y_i, 4)*14.0,
           gsl_vector_get(Y_i, 5)*15.0,
           gsl_vector_get(Y_i, 6)*15.0,
           gsl_vector_get(Y_i, 7)*16.0,
           gsl_vector_get(Y_i, 8)*17.0,
           gsl_vector_get(Y_i, 9)*17.0,
           gsl_vector_get(Y_i,10)*18.0,
           gsl_vector_get(Y_i,11)*18.0,
           gsl_vector_get(Y_i,12)*1.0);

	   // move forward in time and repeat
           t_now += h;
  } while (t_now <= t_stop);

  // de-allocate variables
  gsl_matrix_free(jac);
  gsl_matrix_free(lhs);
  gsl_vector_free(Y_i);
  gsl_vector_free(delta);
  gsl_vector_free(ode_rhs);
  gsl_permutation_free(p);
  return 0;
}
