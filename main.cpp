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
  const double t_stop = 22.0;

  gsl_matrix *jac = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_matrix *lhs = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_vector *Y_i = gsl_vector_calloc(N_ISO);
  gsl_vector *ode_rhs = gsl_vector_calloc(N_ISO);
  gsl_vector *delta = gsl_vector_calloc(N_ISO);
  gsl_permutation *p = gsl_permutation_alloc(N_ISO);

  // mostly H1, with a tiny bit of C12, otherwise the rxn won't go!
  gsl_vector_set(Y_i, N_ISO-1, 0.99);
  gsl_vector_set(Y_i, 1, 8.3e-4);

  printf("%12s %12s\n", "t_now", "Y(H)");
  do {
    error = create_jacobian_matrix(jac, Y_i, T);
    //printf("\nprinting jacobian matrix...\n");
    //gsl_matrix_fprintf(stdout, jac, "%12.3e");

    gsl_matrix_set_identity(lhs);
    //printf("\n\nprinting identity matrix...\n");
    //gsl_matrix_fprintf(stdout, lhs, "%12.3e");
    for (unsigned int i = 0; i < N_ISO; i++) {
      gsl_matrix_set(lhs, i, i, gsl_matrix_get(lhs, i, i) / h);
    }
    error = gsl_matrix_sub(lhs, jac);
    //printf("\n\nprinting 1/h - J...\n");
    //gsl_matrix_fprintf(stdout, lhs, "%12.6e");

    error = build_ode_rhs(ode_rhs, Y_i, T);

    gsl_linalg_LU_decomp(lhs, p, &s);
    gsl_linalg_LU_solve(lhs, p, ode_rhs, delta);

    //gsl_vector_fprintf(stdout, delta, "%12.6f");

    gsl_vector_add(Y_i, delta);
    printf("%12.9f %12.9e\n", t_now, gsl_vector_get(Y_i, N_ISO-1));
    //gsl_vector_fprintf(stdout, Y_i, "%12.3e");
    t_now += h;
  } while (t_now <= t_stop);

  gsl_matrix_free(jac);
  gsl_matrix_free(lhs);
  gsl_vector_free(Y_i);
  gsl_vector_free(delta);
  gsl_vector_free(ode_rhs);
  gsl_permutation_free(p);
  return 0;
}
