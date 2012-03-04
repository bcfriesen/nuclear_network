/* The subscripts on these rate coefficients are hard-coded for the isotopes
 * in the CNO cycle, as defined in my notes. This is a terrible way to label
 * them if one plans to add or remove isotopes from the network in the future,
 * but it's faster to write code this way which is helpful when a deadline
 * looms, as is the case here. */

// i, j = starting products. T = temperature (K)
double lambda_ij (int i, int j, double T) {
  double lambda_ij;
  if (i == 7 && j == 13) {
    lambda_ij = lambda_N15_P_G_O16(T);
  }
  return lambda_ij;
}

/* same, except an extra argument to distinguish a few reactions which have
 * the same starting products but different final products */
double lambda_ij (int i, int j, double T, bool product) {
}

double lambda_N15_P_G_O16 (double T) {
  double T9 = T / 1.0e+09;
  double lambda = 9.78e+08 / pow(T9, 2.0/3.0)
                * exp((-15.251 / pow(T9, 1.0/3.0))
                - pow(T9 / 0.450, 2.0)) * (1.0 + 0.027 * pow(T9, 1.0/3.0)
		+ 0.219 * pow(T9, 2.0/3.0) + 0.042 * T9 + 6.83 * pow(T9, 4.0/3.0)
		+ 3.32 * pow(T9, 5.0/3.0))
		+ 1.11e+04 / pow(T9, 3.0/2.0) * exp(-3.328 / T9)
		+ 1.49e+04 / pow(T9, 3.0/2.0) * exp(-4.665 / T9)
		+ 3.80e+06 / pow(T9, 3.0/2.0) * exp(-11.048 / T9);
}
