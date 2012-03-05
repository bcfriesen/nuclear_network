#include "rate_coeffs.hpp"

#define N_NONZERO_EL 47
// Create sparse Jacobian matrix. For this CNO network the matrix is 72% sparse.
void create_jacobian_matrix (double *, double);

void create_jacobian_matrix (double *Y_i, double T) {
  /* 47 elements are non-zero for this reaction network. These results took
   * about an hour to calculate by hand. */
  
  /* In my notes I started matrix elements at 1, not 0. So these values are
   * universally lower by 1. */
  double mtx_el[N_NONZERO_EL];

  // J_{0, 6}
  mtx_el[0] = lambda_ij(6,12,T,'a')*Y_i[12];
  // J_{0, 9}
  mtx_el[1] = lambda_ij(9,12,T,'a')*Y_i[12];
  // J_{0, 11}
  mtx_el[2] = lambda_ij(11,12,T)*Y_i[12];
  // J_{0, 12}
  mtx_el[3] = lambda_ij(9,12,T,'a')*Y_i[6] + lambda_ij(9,12,T,'a')*Y_i[9]
            + lambda_ij(11,12,T)*Y_i[11];

  // J_{1, 1}
  mtx_el[4] = -lambda_ij(1,12,T)*Y_i[12];
  // J_{1, 6}
  mtx_el[5] = lambda_ij(6,12,T)*Y_i[12];
  // J_{1, 12}
  mtx_el[6] = -lambda_ij(1,12,T)*Y_i[1] + lambda_ij(6,12,T)*Y_i[6];

  // J_{2, 1}
  mtx_el[7] = lambda_ij(1,12,T)*Y_i[12];
  // J_{2, 2}
  mtx_el[8] = -lambda_ij(2);
  // J_{2, 12}
  mtx_el[9] = lambda_ij(1,12,T)*Y_i[1];

  // J_{3, 2}
  mtx_el[10] = lambda_ij(2);
  // J_{3, 3}
  mtx_el[11] = -lambda_ij(3,12,T)*Y_i[12];
  // J_{3, 12}
  mtx_el[12] = -lambda_ij(3,12,T)*Y_i[3];

  // J_{4, 3}
  mtx_el[13] = lambda_ij(3,12,T)*Y_i[12];
  // J_{4, 4}
  mtx_el[14] = -lambda_ij(4,12,T)*Y_i[12];
  // J_{4, 9}
  mtx_el[15] = lambda_ij(9,12,T)*Y_i[12];
  // J_{4, 12}
  mtx_el[16] = -lambda_ij(4,12,T)*Y_i[4] + lambda_ij(3,12,T)*Y_i[3]
             +  lambda_ij(9,12,T)*Y_i[9];

  // J_{5, 4}
  mtx_el[17] = lambda_ij(4,12,T)*Y_i[12];
  // J_{5, 5}
  mtx_el[18] = -lambda_ij(5);
  // J_{5, 12}
  mtx_el[19] = lambda_ij(4,12,T)*Y_i[4];

  // J_{6, 5}
  mtx_el[20] = lambda_ij(5);
  // J_{6, 6}
  mtx_el[21] = -lambda_ij(6,12,T,'a')*Y_i[12] - lambda_ij(6,12,T,'g')*Y_i[12];
  // J_{6, 11}
  mtx_el[22] = lambda_ij(11,12,T)*Y_i[12];
  // J_{6, 12}
  mtx_el[23] = -lambda_ij(6,12,T,'a')*Y_i[6] - lambda_ij(6,12,T,'g')*Y_i[6]
             +  lambda_ij(11,12,T)*Y_i[12];

  // J_{7, 6}
  mtx_el[24] = lambda_ij(6,12,T)*Y_i[12];
  // J_{7, 7}
  mtx_el[25] = -lambda_ij(7,12,T)*Y_i[12];
  // J_{7, 12}
  mtx_el[26] = lambda_ij(6,12,T)*Y_i[6] - lambda_ij(7,12,T)*Y_i[7];

  // J_{8, 7}
  mtx_el[27] = lambda_ij(7,12,T)*Y_i[12];
  // J_{8, 8}
  mtx_el[28] = -lambda_ij(8);
  // J_{8, 12}
  mtx_el[29] = lambda_ij(7,12,T)*Y_i[7];

  // J_{9, 8}
  mtx_el[30] = lambda_ij(8);
  // J_{9, 9}
  mtx_el[31] = -lambda_ij(9,12,T,'g')*Y_i[12] - lambda_ij(9,12,T,'a')*Y_i[12];
  // J_{9, 12}
  mtx_el[32] = -lambda_ij(9,12,T,'g')*Y_i[9] - lambda_ij(9,12,T,'a')*Y_i[9];

  // J_{10, 9}
  mtx_el[33] = lambda_ij(9,12,T,'g')*Y_i[12];
  // J_{10, 10}
  mtx_el[34] = -lambda_ij(10);
  // J_{10, 12}
  mtx_el[35] = lambda_ij(9,12,T,'g')*Y_i[9];

  // J_{11, 10}
  mtx_el[36] = lambda_ij(10);
  // J_{11, 11}
  mtx_el[37] = -lambda_ij(11,12,T)*Y_i[12];
  // J_{11, 12}
  mtx_el[38] = -lambda_ij(11,12,T)*Y_i[11];
  
  // J_{12, 1}
  mtx_el[39] = -lambda_ij(2,12,T)*Y_i[12];
  // J_{12, 3}
  mtx_el[40] = -lambda_ij(3,12,T)*Y_i[12];
  // J_{12, 4}
  mtx_el[41] = -lambda_ij(4,12,T)*Y_i[12];
  // J_{12, 6}
  mtx_el[42] = -lambda_ij(6,12,T,'a')*Y_i[12] - lambda_ij(6,12,T,'g')*Y_i[12];
  // J_{12, 7}
  mtx_el[43] = -lambda_ij(7,12,T)*Y_i[12];
  // J_{12, 9}
  mtx_el[44] = -lambda_ij(9,12,T,'a')*Y_i[12] - lambda_ij(9,12,T,'g')*Y_i[12];
  // J_{12, 11}
  mtx_el[45] = -lambda_ij(11,12,T)*Y_i[12];
  // J_{12, 12}
  mtx_el[46] = -lambda_ij(1,12,T)*Y_i[1] - lambda_ij(3,12,T)*Y_i[3]
             -  lambda_ij(4,12,T)*Y_i[4] - lambda_ij(6,12,T,'a')*Y_i[6]
	     -  lambda_ij(6,12,T,'g')*Y_i[6] - lambda_ij(7,12,T)*Y_i[7]
	     -  lambda_ij(9,12,T,'g')*Y_i[9] - lambda_ij(9,12,T,'a')*Y_i[9]
	     -  lambda_ij(11,12,T)*Y_i[11];
  return;
}
