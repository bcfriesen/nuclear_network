#!/bin/bash

icpc -I /home5/baker/friesen/gsl-1.15/intel-12.1.2/include \
     -L /home5/baker/friesen/gsl-1.15/intel-12.1.2/lib \
     -lgsl -lgslcblas \
     -mkl -O3 \
     -o run \
     rate_coeffs.cpp jacobian.cpp ode_rhs.cpp main.cpp

#icpc -I /home5/baker/friesen/gsl-1.15/intel-12.1.2/include \
#     -L /home5/baker/friesen/gsl-1.15/intel-12.1.2/lib \
#     -lgsl -lgslcblas \
#     -mkl -O3 \
#     -o test2 rate_coeffs.cpp \
#     test2.cpp
