#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "jacobian.hpp"
#include "ode_rhs.hpp"
#include "rate_coeffs.hpp"
#include "global.hpp"

int main() {
  // temperature (constant throughout). units: K
  double T = 15.0e+06;
  // density (constant throughout). units: g/cm^3
  double rho = 150.0;
  // units: g/mol
  const double molar_mass[nvar] = {4.002602, 12.0, 13.005738609, 13.00335483778,
  14.00307400478, 15.003065617, 15.00010889823, 15.99491461956, 17.002095237,
  16.999131703, 18.000937956, 17.999161001, 1.00794};
  /* initial time step (sec). This is just an initial guess. The
   * time-stepper will fix it when it starts integrating. */
  double h = 1.0e0;
  // initial and final times. units: sec
  double t_now = 0.0, t_stop = 1.0e+22;
  // absolute and relative error requirements
  const double eps_abs = 1.0e-8, eps_rel = 0.0;

  // number abundances of isotopes. units: mol/cm^3
  double y[nvar];
  /* jacobian matrix. the integrator needs this. fortunately it's
   * analytic so calculating it is very fast */
  double dfdy[nvar][nvar];
  /* time derivatives of the RHS. the integrator needs this. there is
   * no explicit time dependence in this system of ODEs so these will
   * all be zero */
  double dfdt[nvar];

  /* initial abundances. these are sort of arbitrary. I assume the
   * environment is the core of a young star, so 99% H1 (by mass) and
   * 1% C12 (we only need a tiny bit of C12 to start the reaction) */
  for (unsigned int i = 0; i < nvar; ++i) {
    /* the integrator doesn't like "true" zeros very much, so we use
     * tiny positive numbers instead. the stuff in parentheses
     * converts mass fraction to mol/cm^3 */
    y[i] = 1.0e-20 * (rho / molar_mass[i]);
  }
  // set H1 and C12 by hand
  y[12] = 0.99 * (rho / molar_mass[12]);
  y[ 1] = 0.01 * (rho / molar_mass[ 1]);

  /* declare integration technology. All this junk is built in to the
   * GNU Scientific Library. I'm using a Bulirsch-Stoer integration
   * method, which is a variable-order (from 7th to 15th) method
   * designed specifically for integrating extremely stiff systems of
   * ODEs, like nuclear networks. the benefit of using a sophisticated
   * algorithm like this is that it adjusts the time steps based on
   * whether something interesting is happening or not.  In "Numerical
   * Recipes" they compare Bulirsch-Stoer to a plain-jane 4th-order
   * Runge-Kutta scheme for integrating some system of ODEs, and the
   * Runge-Kutta method takes something like 50,000 time steps to
   * solve the equations, whereas B-S took only 29. */
  const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_bsimp;
  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(step_type, nvar);
  gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(nvar);
  gsl_odeiv2_system sys = {ode_rhs, jacobian, nvar, &T};

  FILE *fp;
  fp = fopen("results.dat", "w");
  // print column headers
  fprintf(fp, "%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s"
         " %15s\n", "tnow", "he4", "c12", "n13", "c13", "n14", "o15",
         "n15", "o16", "f17", "o17", "f18", "o18", "h1");
  // continue loop until we reach t_stop
  while (t_now < t_stop) {
    // the integrator will update t_now after each time step.
    int status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys,
                 &t_now, t_stop, &h, y);
    // quit if there's an error
    if (status != GSL_SUCCESS) break;
    // print isotope mass fractions at each time step
    fprintf(fp, "%15.4e %15.4e %15.4e %15.4e %15.4e %15.4e %15.4e %15.4e %15.4e"
           " %15.4e %15.4e %15.4e %15.4e %15.4e\n",
	   t_now,
	   y[0] / (rho / molar_mass[0]),
	   y[1] / (rho / molar_mass[1]),
	   y[2] / (rho / molar_mass[2]),
	   y[3] / (rho / molar_mass[3]),
	   y[4] / (rho / molar_mass[4]),
	   y[5] / (rho / molar_mass[5]),
	   y[6] / (rho / molar_mass[6]),
	   y[7] / (rho / molar_mass[7]),
	   y[8] / (rho / molar_mass[8]),
	   y[9] / (rho / molar_mass[9]),
	   y[10] / (rho / molar_mass[10]),
	   y[11] / (rho / molar_mass[11]),
	   y[12] / (rho / molar_mass[12]));
  }

  gsl_odeiv2_step_free(step);
  gsl_odeiv2_control_free(control);
  gsl_odeiv2_evolve_free(evolve);
  fclose(fp);
  return 0;
}
