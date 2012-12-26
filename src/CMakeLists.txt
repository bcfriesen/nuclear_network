SET (nuclear_network_SOURCES
jacobian.c
main.c
ode_rhs.c
rate_coeffs.c
)

FIND_PACKAGE(GSL)
IF (GSL_FOUND)
  SET (CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${GSL_C_FLAGS})
ENDIF (GSL_FOUND)

ADD_EXECUTABLE (nuclear_network ${nuclear_network_SOURCES})
TARGET_LINK_LIBRARIES(nuclear_network
    gsl
    gslcblas
    )
