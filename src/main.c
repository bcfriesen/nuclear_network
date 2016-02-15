#define MAIN_FILE

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "ode_rhs.h"
#include "rate_coeffs.h"
#include "jacobian.h"
#include "param.h"

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

const int n_iso = 13;           // # of isotopes to evolve

// # of grid points in each direction
const int nx = 4;
const int ny = 4;
const int nz = 4;

/* molar masses of each isotope. used to convert from mass fraction to molar
 * number abundance. units: g/mol */
double *molar_mass;

unsigned int i, j;

// state data stored in each grid cell
struct grid_cell
{
  struct param params;
  double *y;                    // molar (number) fraction of isotopes
  double *dfdy;                 // Jacobian
  /* time derivatives of the RHS. the integrator needs this. there is
   * no explicit time dependence in this system of ODEs so these will
   * all be zero */
  double *dfdt;
};


void
initial_conditions (struct grid_cell *cell)
{
  // temperature (constant throughout). units: K
  cell->params.T = 25.0e+06;
  // density (constant throughout). units: g/cm^3
  cell->params.rho = 150.0;
  // number of isotopes to include in network
  cell->params.n_iso = 13;

  /* set initial abundances. these are sort of arbitrary. I assume the
   * environment is the core of a young star, so 99% H1 (by mass) and
   * 1% C12 (we only need a tiny bit of C12 to start the reaction) */
  for (j = 0; j < cell->params.n_iso; ++j)
    {
      /* the integrator doesn't like "true" zeros very much, so we use
       * tiny positive numbers instead. the stuff in parentheses
       * converts mass fraction to mol/cm^3 */
      cell->y[j] = 1.0e-20 * (cell->params.rho / molar_mass[j]);
    }
  // set H1 and C12 by hand
  cell->y[12] = 0.99 * (cell->params.rho / molar_mass[12]);
  cell->y[1] = 0.01 * (cell->params.rho / molar_mass[1]);
}


int
main ()
{

  /* initial time step (sec). This is just an initial guess. The
   * time-stepper will fix it when it starts integrating. */
  double h = 1.0e-8;
  /* initial and final times. units: sec. the abundances for this problem
   * should evolve on stellar evolution timescales. for reference,
   * 1 Gyr ~ 3e16 sec */
  double t_now = 0.0, t_stop = 1.0e+10;
  /* absolute and relative error requirements for the integrator. smaller means
   * better precision but more computation time */
  const double eps_abs = 1.0e-8, eps_rel = 0.0;

  // total # of grid points
  int tot_grid_pts = nx * ny * nz;

  molar_mass = malloc (n_iso * sizeof (double));
  molar_mass[0] = 4.002602;     // He4
  molar_mass[1] = 12.0;         // C12
  molar_mass[2] = 13.005738609; // N13
  molar_mass[3] = 13.00335483778;       // C13
  molar_mass[4] = 14.00307400478;       // N14
  molar_mass[5] = 15.003065617; // O15
  molar_mass[6] = 15.00010889823;       // N15
  molar_mass[7] = 15.99491461956;       // O16
  molar_mass[8] = 17.002095237; // F17
  molar_mass[9] = 16.999131703; // O17
  molar_mass[10] = 18.000937956;        // F18
  molar_mass[11] = 17.999161001;        // O18
  molar_mass[12] = 1.00794;     // H1

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
  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc (step_type, n_iso);
  // set absolute and relative error tolerances
  gsl_odeiv2_control *control = gsl_odeiv2_control_y_new (eps_abs, eps_rel);
  // set number of ODEs to solve
  gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc (n_iso);
  /* the integrator needs to know the RHS of the ODEs (ode_rhs), the
   * Jacobian matrix (jacobian), the number of ODEs it's going to
   * solve (n_iso), and any additional parameters (just temperature in
   * this case) */

  unsigned int i, j;

  // create grid
  struct grid_cell grid[tot_grid_pts];
  for (i = 0; i < tot_grid_pts; ++i)
    {
      grid[i].dfdy = malloc (n_iso * n_iso * sizeof (double));
      grid[i].dfdt = malloc (n_iso * sizeof (double));
      grid[i].y = malloc (n_iso * sizeof (double));
    }

  // Set up initial conditions.
  for (i = 0; i < tot_grid_pts; ++i)
    {
      initial_conditions (&(grid[i]));
    }

  /* Duplicates of the current time step and step size. The ODE driver updates
   * these in-place but we want to use the same values to be used on entry at
   * every grid cell. If we don't keep "pristine" copies of the original time
   * step and step size prior to each grid sweep, then the input value to the
   * ODE driver at each cell will be the output from the integration in the
   * previous cell, which we don't want. */
  double t_now_tmp, h_tmp;

  while (t_now < t_stop)
    {
      printf ("t = %e\n", t_now);

      // Grid sweep.
      for (i = 0; i < tot_grid_pts; ++i)
        {

          gsl_odeiv2_system sys =
            { ode_rhs, jacobian, n_iso, &(grid[i].params) };

          /* Integrate the equations at time t_now and take a step forward. The
           * integrator will update the current time in-place, but we don't save that
           * value until after the last grid point is done. */

          t_now_tmp = t_now;
          h_tmp = h;
          int status = gsl_odeiv2_evolve_apply (evolve, control, step, &sys,
                                                &t_now_tmp, t_stop, &h_tmp,
                                                grid[i].y);
          if (status != GSL_SUCCESS)
            break;

          /* Kill an isotope if its mass fraction drops below some really
           * small value. This helps the integrator move a little faster
           * because otherwise it tries to resolve changes at like 1.0e-58,
           * which grids the whole process to a halt. */
          for (j = 0; j < grid[i].params.n_iso; ++j)
            {
              if (grid[i].y[j] / (grid[i].params.rho / molar_mass[j]) <
                  1.0e-20)
                {
                  grid[i].y[j] = 0.0;
                }
            }

        }

      /* Now that the grid sweep is done, save the last value of the new current
       * time and time step. The results from every cell in this example should be
       * bitwise identical so we can just save the very last one. */
      t_now = t_now_tmp;
      h = h_tmp;

    }

  // free GSL memory
  gsl_odeiv2_step_free (step);
  gsl_odeiv2_control_free (control);
  gsl_odeiv2_evolve_free (evolve);

  // free grid data
  for (i = 0; i < tot_grid_pts; ++i)
    {
      free (grid[i].y);
      free (grid[i].dfdy);
      free (grid[i].dfdt);
    }

  free (molar_mass);

  return 0;
}
