#include "rate_coeffs.hpp"

int ode_rhs(double t, const double y[], double dydt[], void *params);
  double temp = *(double *)params {
  dydt[ 0] =  lambda_ij( 6, 12, temp, 'a') * y[ 6] * y[12]
            + lambda_ij( 9, 12, temp, 'a') * y[ 9] * y[12]
            + lambda_ij(11, 12, temp     ) * y[11] * y[12];
  dydt[ 1] = -lambda_ij( 1, 12, temp     ) * y[ 1] * y[12]
            + lambda_ij( 6, 12, temp, 'a') * y[ 6] * y[12];
  dydt[ 2] =  lambda_ij( 1, 12, temp     ) * y[ 1] * y[12]
            - lambda_ij( 2               ) * y[ 2];
  dydt[ 3] =  lambda_ij( 2               ) * y[ 2]
            - lambda_ij( 3, 12, temp     ) * y[ 3] * y[12];
  dydt[ 4] = -lambda_ij( 4, 12, temp     ) * y[ 4] * y[12]
            + lambda_ij( 3, 12, temp     ) * y[ 3] * y[12]
	    + lambda_ij( 9, 12, temp, 'a') * y[ 9] * y[12];
  dydt[ 5] = -lambda_ij( 5               ) * y[ 5]
            + lambda_ij( 4, 12, temp     ) * y[ 4] * y[12];
  dydt[ 6] =  lambda_ij( 5               ) * y[ 5]
            - lambda_ij( 6, 12, temp, 'a') * y[ 6] * y[12]
	    - lambda_ij( 6, 12, temp, 'g') * y[ 6] * y[12]
	    + lambda_ij(11, 12, temp     ) * y[11] * y[12];
  dydt[ 7] =  lambda_ij( 6, 12, temp, 'g') * y[ 6] * y[12]
            - lambda_ij( 7, 12, temp     ) * y[ 7] * y[12];
  dydt[ 8] =  lambda_ij( 7, 12, temp     ) * y[ 7] * y[12]
            - lambda_ij( 8               ) * y[ 8];
  dydt[ 9] =  lambda_ij( 8               ) * y[ 8]
            - lambda_ij( 9, 12, temp, 'g') * y[ 9] * y[12]
	    - lambda_ij( 9, 12, temp, 'a') * y[ 9] * y[12];
  dydt[10] =  lambda_ij( 9, 12, temp, 'g') * y[ 9] * y[12]
            - lambda_ij(10               ) * y[10];
  dydt[11] =  lambda_ij(10               ) * y[10]
            - lambda_ij(11, 12, temp     ) * y[11] * y[12];
  dydt[12] = -lambda_ij( 1, 12, temp     ) * y[ 1] * y[12]
            - lambda_ij( 3, 12, temp     ) * y[ 3] * y[12]
	    - lambda_ij( 4, 12, temp     ) * y[ 4] * y[12]
	    - lambda_ij( 6, 12, temp, 'a') * y[ 6] * y[12]
	    - lambda_ij( 6, 12, temp, 'g') * y[ 6] * y[12]
	    - lambda_ij( 7, 12, temp     ) * y[ 7] * y[12]
	    - lambda_ij( 9, 12, temp, 'a') * y[ 9] * y[12]
	    - lambda_ij( 9, 12, temp, 'g') * y[ 9] * y[12]
	    - lambda_ij(11, 12, temp     ) * y[11] * y[12];
  return 0;
}
