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

  gsl_matrix *jac = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_matrix *lhs = gsl_matrix_calloc(N_ISO, N_ISO);
  gsl_vector *Y_i = gsl_vector_calloc(N_ISO);
  gsl_vector *ode_rhs = gsl_vector_calloc(N_ISO);
  gsl_vector *delta = gsl_vector_calloc(N_ISO);
  gsl_permutation *p = gsl_permutation_alloc(N_ISO);

  // mostly H1, with a tiny bit of C12, otherwise the rxn won't go!
  gsl_vector_set(Y_i, N_ISO-1, 0.99);
  gsl_vector_set(Y_i, 1, 8.3e-4);

  printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
  "t_now", "Y(He4)", "Y(C12)", "Y(N13)", "Y(C13)",
  "Y(N14)", "Y(O15)", "Y(N15)", "Y(O16)", "Y(F17)", "Y(O17)", "Y(F18)",
  "Y(O18)", "Y(H1)");
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
    printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",
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
