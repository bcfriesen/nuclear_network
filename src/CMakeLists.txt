SET (nuclear_network_SOURCES
jacobian.c
main.c
ode_rhs.c
rate_coeffs.c
)

ADD_EXECUTABLE (nuclear_network ${nuclear_network_SOURCES})
TARGET_LINK_LIBRARIES(nuclear_network
    gsl
    gslcblas
    m
    )
