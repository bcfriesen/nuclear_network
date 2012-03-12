#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "jacobian.hpp"
#include "ode_rhs.hpp"
#include "rate_coeffs.hpp"

size_t nvar = 13;

int main() {
  const double T = 15.0e+06;

  unsigned int i;
  double temp = 15.0e+06;

  double dfdy[nvar][nvar], y[nvar];

  /*
  printf("%12s %12s\n", "T (K)", "lambda");
  for (i = 0; i < 100; i++) {
    temp = 1.0e+06 + (double)i * (1.0e+07 - 1.0e+06) / 100.0;
    printf("%12.3f %12.3e\n", temp, lambda_ij(6, 12, temp, 'g'));
  }
  */

  int s, error;

  // time step (sec)
  double h = 1.0e-6;
  double t_now = 0.0, t_stop = 22.0;

  const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_bsimp;
  gsl_odeiv2_step *step       = gsl_odeiv2_step_alloc(step_type, nvar);
  gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1.0e-6, 0.0);
  gsl_odeiv2_evolve *evolve   = gsl_odeiv2_evolve_alloc(nvar);

  gsl_odeiv2_system sys = {ode_rhs, jacobian, nvar, &temp};

  // pure H
  y[12] = 0.99;
  y[ 1] = 8.3e-4;

  printf("%12s %12s\n", "t_now", "Y(H)");
  while (t_now < t_stop) {
    int status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys,
                 &t_now, t_stop, &h, y);
    if (status != GSL_SUCCESS) break;
    printf("%12.4e %12.4e\n", t_now, y[0]);
  }

  gsl_odeiv2_evolve_free(evolve);
  gsl_odeiv2_control_free(control);
  gsl_odeiv2_step_free(step);
  return 0;
}
