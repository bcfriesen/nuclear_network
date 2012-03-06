#include <math.h>
#include <gsl/gsl_const_num.h>
#include "rate_coeffs.hpp"

/* The subscripts on these rate coefficients are hard-coded for the isotopes
 * in the CNO cycle, as defined in my notes. This is a terrible way to label
 * them if one plans to add or remove isotopes from the network in the future,
 * but it's faster to write code this way which is helpful when a deadline
 * looms, as is the case here. */

// All 2-body rates taken from Caughlan & Fowler (1988)
//
// All beta-decay rates taken from Wolfram|Alpha. Remember that these rates
// are temperature-independent.

// i, j = starting products. T = temperature (K)

/* In my notes I started all matrix/vector elements at 1, not 0, so the values
 * here will be universally lower by 1. */
double lambda_ij (int i, int j, double T) {
  double lambda_ij;
  if        (i == 11 && j == 12) {
    lambda_ij = lambda_O18_P_A_N15(T);
  } else if (i == 1 && j == 12) {
    lambda_ij = lambda_C12_P_G_N13(T);
  } else if (i == 3 && j == 12) {
    lambda_ij = lambda_C13_P_G_N14(T);
  } else if (i == 4 && j == 12) {
    lambda_ij = lambda_N14_P_G_O15(T);
  } else if (i == 7 && j == 12) {
    lambda_ij = lambda_O16_P_G_F17(T);
  } else {
    lambda_ij = 0.0;
  }
  return lambda_ij * GSL_CONST_NUM_AVOGADRO;
}

/* same, except an extra argument to distinguish a few reactions which have
 * the same starting products but different final products */
double lambda_ij (int i, int j, double T, char product) {
  double lambda_ij;
  if        (i == 6 && j == 12 && product == 'a') {
    lambda_ij = lambda_N15_P_A_C12(T);
  } else if (i == 6 && j == 12 && product == 'g') {
    lambda_ij = lambda_N15_P_G_O16(T);
  } else if (i == 9 && j == 12 && product == 'a') {
    lambda_ij = lambda_O17_P_A_N14(T);
  } else if (i == 9 && j == 12 && product == 'g') {
    lambda_ij = lambda_O17_P_G_F18(T);
  } else {
    lambda_ij = 0.0;
  }
  return lambda_ij * GSL_CONST_NUM_AVOGADRO;
}

// beta-decay reactions aren't 2 body and so only require 1 argument
double lambda_ij (int i) {
  double lambda_ij;
  if        (i == 2) {
    lambda_ij = lambda_N13_e_nu();
  } else if (i == 5) {
    lambda_ij = lambda_O15_e_nu();
  } else if (i == 8) {
    lambda_ij = lambda_F17_e_nu();
  } else if (i == 10) {
    lambda_ij = lambda_F18_e_nu();
  } else {
    lambda_ij = 0.0;
  }
  return lambda_ij * GSL_CONST_NUM_AVOGADRO;
}

// lambda_{6, 12, alpha}
double lambda_N15_P_A_C12 (double T) {
  // TODO: figure what this stupid fudge factor is
  double T9 = T / 1.0e+09;
  double T913 = pow(T9, 1.0/3.0), T923 = pow(T9, 2.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0), T912 = pow(T9, 1.0/2.0);
  // mmmm, double fudge sounds delicious...
  double fudge = 0.5;
  double lambda = 1.08e+12 / T923 * exp(-15.251 / T913 - pow(T9 / 0.522, 2.0))
  * (1.0 + 0.027 * T913 + 2.62 * T923 + 0.501 * T9 + 5.36 * T943 + 2.60 * T953)
  + 1.19e+08 / T932 * exp(-3.676 / T9) + 5.41e+08 / T912 * exp(-8.926 / T9)
  + fudge * 4.72e+08 / T932 * exp(-7.721 / T9) + 2.20e+09 / T932
  * exp(-11.418 / T9);
  return lambda;
}

// lambda_{9, 12, alpha}
double lambda_O17_P_A_N14 (double T) {
  double T9 = T / 1.0e+09;
  // TODO: check out the weird message in CF88 about fudge2
  double fudge = 0.5, fudge2 = 0.5;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0);
  double lambda = 1.53e+07 / T923 * exp(-16.712 / T913
  - pow(T9 / 0.565, 2.0)) * (1.0 + 0.025 * T913 + 5.39 * T923 + 0.940 * T9
  + 13.5 * T943 + 5.98 * T953) + fudge * (4.81e+10 * T9 * exp(-16.712 / T913
  - pow(T9 / 0.040, 2.0)) + 5.05e-05 / T932 * exp(-0.723 / T9)) + fudge2
  * 1.31e+01 / T932 * exp(-1.961 / T9);
  return lambda;
}

// lambda_{11, 12}
double lambda_O18_P_A_N15 (double T) {
  double T9 = T / 1.0e+09;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0);
  double lambda = 3.63e+11 / T923 * exp(-16.729 / T913 - pow(T9 / 1.361, 2.0))
  * (1.0 + 0.025 * T913 + 1.88 * T923 + 0.327 * T9 + 4.66 * T943 + 2.06 * T953)
  + 9.90e-14 / T932 * exp(-0.231 / T9) + 2.66e+04 / T932 * exp(-1.670 / T9)
  + 2.41e+09 / T932 * exp(-7.638 / T9) + 1.46e+09 / T9 * exp(-8.310 / T9);
  return lambda;
}

// lambda_{1, 12}
double lambda_C12_P_G_N13 (double T) {
  double T9 = T / 1.0e+09;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0);
  double lambda = 2.04e+07 / T923 * exp(-13.690 / T913 - pow(T9 / 1.500, 2.0))
  * (1.0 + 0.030 * T913 + 1.19 * T923 + 0.254 * T9 + 2.06 * T943 + 1.12 * T953)
  + 1.08e+05 / T932 * exp(-4.925 / T9) + 2.15e+05 / T932 * exp(-18.179 / T9);
  return lambda;
}

// lambda_{3, 12}
double lambda_C13_P_G_N14 (double T) {
  double T9 = T / 1.0e+09;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0), T965 = pow(T9, 6.0/5.0);
  double lambda = 8.01e+07 / T923 * exp(-13.717 / T913 - pow(T9 / 2.000, 2.0))
  * (1.0 + 0.030 * T913 + 0.958 * T923 + 0.204 * T9 + 1.39 * T943 + 0.753 * T953)
  + 1.21e+06 / T965 * exp(-5.701 / T9);
  return lambda;
}

// lambda_{4, 12}
double lambda_N14_P_G_O15 (double T) {
  double T9 = T / 1.0e+09;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0), T965 = pow(T9, 6.0/5.0);
  double lambda = 4.90e+07 / T923 * exp(-15.202 / T913 - pow(T9 / 1.191, 2.0))
  * (1.0 + 0.027 * T913 - 0.778 * T923 - 0.149 * T9 + 0.261 * T943
  * + 0.127 * T953) + 2.37e+03 / T932 * exp(-3.011 / T9)
  + 2.19e+04 * exp(-12.530/ T9);
  return lambda;
}

// lambda_{6, 12, gamma}
double lambda_N15_P_G_O16 (double T) {
  double T9 = T / 1.0e+09;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0), T965 = pow(T9, 6.0/5.0);
  double lambda = 9.78e+08 / pow(T9, 2.0/3.0) * exp((-15.251 / pow(T9, 1.0/3.0))
  - pow(T9 / 0.450, 2.0)) * (1.0 + 0.027 * pow(T9, 1.0/3.0)
  + 0.219 * pow(T9, 2.0/3.0) + 0.042 * T9 + 6.83 * pow(T9, 4.0/3.0)
  + 3.32 * pow(T9, 5.0/3.0)) + 1.11e+04 / pow(T9, 3.0/2.0) * exp(-3.328 / T9)
  + 1.49e+04 / pow(T9, 3.0/2.0) * exp(-4.665 / T9)
  + 3.80e+06 / pow(T9, 3.0/2.0) * exp(-11.048 / T9);
  return lambda;
}

// lambda_{7, 12}
double lambda_O16_P_G_F17 (double T) {
  double T9 = T / 1.0e+09;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0), T965 = pow(T9, 6.0/5.0);
  double lambda = 1.50e+08 / (T923 * (1.0 + 2.13 * (1.0 - exp(-0.728 * T923))))
  * exp(-16.692 / T913);
  return lambda;
}

// lambda_{9, 12, gamma}
double lambda_O17_P_G_F18 (double T) {
  double T9 = T / 1.0e+09;
  double T923 = pow(T9, 2.0/3.0), T913 = pow(T9, 1.0/3.0),
         T943 = pow(T9, 4.0/3.0), T953 = pow(T9, 5.0/3.0),
	 T932 = pow(T9, 3.0/2.0), T965 = pow(T9, 6.0/5.0);
  double T9A = T9 / (1.0 + 2.69 * T9);
  double T9A56 = pow(T9A, 5.0/6.0), T9A13 = pow(T9A, 1.0/3.0);
  // TODO: figure out this fudge factor
  double fudge = 0.5;
  double lambda = 7.97e+07 * T9A56 / T932 * exp(-16.712 / T9A13)
  + 1.51e+08 / T923 * exp(-16.712 / T913) * (1.0 + 0.025 * T913 - 0.051 * T923
  - 8.82e-03 * T9) + 1.56e+05 / T9 * exp(-6.272 / T9) + fudge * 1.31e+01 / T932
  * exp(-1.961 / T9);
  return lambda;
}

// lambda_{2 beta+}
double lambda_N13_e_nu () {
  // half-life (sec)
  double t12 = 9.965 * 60.0;
  /* I'm not sure how to deal with the time dependence of these beta-decay
   * rates, so I set delta_t to 1 sec. I don't know if that's right but these
   * decays are practically instantaneous compared to the 2-body reactions so
   * hopefully it won't screw everything up... */
  double delta_t = 1.0;
  /* Spontaneous decay rates have the form
   *       N(t) = N_0 * (1/2) ^ (t/t12)
   * where N_0 is the starting number density, t is the elapsed time and t12
   * is the half-life. I wanted to write this in the same way I wrote all the
   * other lambda's, that is,
   *       N = N_0 \lambda_decay
   * No idea if this is right though, particularly because this lambda doesn't
   * have units of a rate. But maybe by setting delta_t = 1 sec we can sort of
   * implicitly add a s^(-1) unit to the lambda? */
  double lambda = pow(1.0/2.0, delta_t/t12);
  return lambda;
}

// lambda_{5 beta+}
double lambda_O15_e_nu () {
  double t12 = 122.24;
  double delta_t = 1.0;
  double lambda = pow(1.0/2.0, delta_t/t12);
  return lambda;
}

// lambda_{8 beta+}
double lambda_F17_e_nu () {
  double t12 = 64.49;
  double delta_t = 1.0;
  double lambda = pow(1.0/2.0, delta_t/t12);
  return lambda;
}

// lambda_{10 beta+}
double lambda_F18_e_nu () {
  double t12 = 109.771 * 60;
  double delta_t = 1.0;
  double lambda = pow(1.0/2.0, delta_t/t12);
  return lambda;
}
