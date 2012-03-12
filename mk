#!/bin/bash

GSL_DIR=/home5/baker/friesen/gsl-1.15/intel-12.1.2

icpc -I ${GSL_DIR}/include \
     -L ${GSL_DIR}/lib -lgsl -lgslcblas \
     -mkl -O3 \
     -o run \
     rate_coeffs.cpp jacobian.cpp ode_rhs.cpp main.cpp
