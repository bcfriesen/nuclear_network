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

  /*
  printf("%12s %12s\n", "T (K)", "lambda");
  for (i = 0; i < 100; i++) {
    temp = 1.0e+06 + (double)i * (1.0e+07 - 1.0e+06) / 100.0;
    printf("%12.3f %12.3e\n", temp, lambda_ij(6, 12, temp, 'g'));
  }
  */

  int s, error;

  // time step (sec)
  double h = 0.01;
  double t_now = 0.0, t_stop = 22.0;

  gsl_matrix *jac = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_matrix *lhs = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_vector *Y_i = gsl_vector_calloc(N_ISO);
  gsl_vector *ode_rhs = gsl_vector_calloc(N_ISO);
  gsl_vector *delta = gsl_vector_calloc(N_ISO);
  gsl_permutation *p = gsl_permutation_alloc(N_ISO);

  // pure H
  gsl_vector_set(Y_i, N_ISO-1, 0.99);
  gsl_vector_set(Y_i, 1, 8.3e-4);

  //printf("%12s %12s\n", "t_now", "Y(H)");
  do {
    error = create_jacobian_matrix(jac, Y_i, T);
    //gsl_matrix_fprintf(stdout, jac, "%12.6f");
    //printf("\n");
    gsl_matrix_set_identity(lhs);
    for (unsigned int i = 0; i < N_ISO; i++) {
      gsl_matrix_set(lhs, i, i, gsl_matrix_get(lhs, i, i) / h);
    }
    error = gsl_matrix_sub(lhs, jac);
    //gsl_matrix_fprintf(stdout, lhs, "%12.6e");

    /*
    printf("before calling build_ode_rhs...\n");
    gsl_vector_fprintf(stdout, ode_rhs, "%12.3e");
    error = build_ode_rhs(ode_rhs, Y_i, T);
    printf("after calling build_ode_rhs...\n");
    gsl_vector_fprintf(stdout, ode_rhs, "%12.3e");
    */

    gsl_linalg_LU_decomp(lhs, p, &s);
    gsl_linalg_LU_solve(lhs, p, ode_rhs, delta);

    //gsl_vector_fprintf(stdout, delta, "%12.6f");

    //gsl_vector_add(Y_i, delta);
    //printf("%12.3f %12.3e\n", t_now, gsl_vector_get(Y_i, N_ISO-1));
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
