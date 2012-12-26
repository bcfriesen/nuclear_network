#define MAIN_FILE

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "ode_rhs.h"
#include "rate_coeffs.h"
#include "global.h"
#include "jacobian.h"

/* A simple CNO nuclear network solver. Thanks Dick Henry for making
 * these projects really open ended! I probably would not have learned
 * nearly as much about nuclear networks otherwise. */

/* Because I'm pretty myopic in general I hard-wired the isotopes into
 * the ODEs, Jacobian, etc. A better way to do this would be to create
 * integer variables called ihe4, ic12, etc. and assign the indices to
 * them. Then if I wanted to add more isotopes later it would be quite
 * simple. Oh well, the deadline for this project is right around the
 * corner. The isotope codes are:
 *
 * He4 =  0
 * C12 =  1
 * N13 =  2
 * C13 =  3
 * N14 =  4
 * O15 =  5
 * N15 =  6
 * O16 =  7
 * F17 =  8
 * O17 =  9
 * F18 = 10
 * O18 = 11
 * H1  = 12
 */

int main() {
  // temperature (constant throughout). units: K
  double T = 25.0e+06;
  // density (constant throughout). units: g/cm^3
  double rho = 150.0;
  printf("%18s %12.4e\n", "TEMPERATURE:", T);
  printf("%18s %12.4e\n", "MASS DENSITY:", rho);
  /* molar masses of each isotope. used to convert from mass fraction to
   * molar number abundance. units: g/mol */
  double molar_mass[nvar];
  /* initial time step (sec). This is just an initial guess. The
   * time-stepper will fix it when it starts integrating. */
  double h = 1.0e-8;
  /* initial and final times. units: sec. the abundances for this problem
   * should evolve on stellar evolution timescales. for reference,
   * 1 Gyr ~ 3e16 sec */
  double t_now = 0.0, t_stop = 1.0e+22;
  /* absolute and relative error requirements for the integrator. smaller means
   * better precision but more computation time */
  const double eps_abs = 1.0e-8, eps_rel = 0.0;

  // number abundances of isotopes. units: mol/cm^3
  double y[nvar];
  /* Jacobian matrix. the integrator needs this. fortunately it's
   * analytic so calculating it is very fast */
  double dfdy[nvar][nvar];
  /* time derivatives of the RHS. the integrator needs this. there is
   * no explicit time dependence in this system of ODEs so these will
   * all be zero */
  double dfdt[nvar];
  // loops
  unsigned int i;

  molar_mass[0] = 4.002602;
  molar_mass[1] = 12.0;
  molar_mass[2] = 13.005738609;
  molar_mass[3] = 13.00335483778;
  molar_mass[4] = 14.00307400478;
  molar_mass[5] = 15.003065617;
  molar_mass[6] = 15.00010889823;
  molar_mass[7] = 15.99491461956;
  molar_mass[8] = 17.002095237;
  molar_mass[9] = 16.999131703;
  molar_mass[10] = 18.000937956;
  molar_mass[11] = 17.999161001;
  molar_mass[12] = 1.00794;

  /* set initial abundances. these are sort of arbitrary. I assume the
   * environment is the core of a young star, so 99% H1 (by mass) and
   * 1% C12 (we only need a tiny bit of C12 to start the reaction) */
  for (i = 0; i < nvar; ++i) {
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
   * method (the "bsimp" in the first line stands for "Bulirsch-Stoer
   * implicit"; other options are "rk4" for 4th-order Runge-Kutta,
   * etc.), which is a variable-order (from 5th to 15th) method
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
  // set absolute and relative error tolerances
  gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  // set number of ODEs to solve
  gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(nvar);
  /* the integrator needs to know the RHS of the ODEs (ode_rhs), the
   * Jacobian matrix (jacobian), the number of ODEs it's going to
   * solve (nvar), and any additional parameters (just temperature in
   * this case) */
  gsl_odeiv2_system sys = {ode_rhs, jacobian, nvar, &T};

  // pointer for writing output to a file
  FILE *fp;
  fp = fopen("results.dat", "w");
  // print column headers
  fprintf(fp, "%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s"
         " %15s\n", "tnow", "he4", "c12", "n13", "c13", "n14", "o15",
         "n15", "o16", "f17", "o17", "f18", "o18", "h1");
  // continue loop until we reach t_stop
  while (t_now < t_stop) {
    /* integrate the equations at time t_now and take a step forward
     * (t_now will be updated automatically) */
    int status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys,
                 &t_now, t_stop, &h, y);
    // quit if there's an error
    if (status != GSL_SUCCESS) break;
    /* Kill an isotope if its mass fraction drops below some really
     * small value. This helps the integrator move a little faster
     * because otherwise it tries to resolve changes at like 1.0e-58,
     * which is pointless. */
    for (i = 0; i < nvar; ++i) {
      if (y[i] / (rho / molar_mass[i]) < 1.0e-20) y[i] = 0.0;
    }
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

  // free pointers
  gsl_odeiv2_step_free(step);
  gsl_odeiv2_control_free(control);
  gsl_odeiv2_evolve_free(evolve);
  // close file
  fclose(fp);
  return 0;
}
