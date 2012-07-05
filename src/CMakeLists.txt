SET (nuclear_network_SOURCES
jacobian.cpp
main.cpp
ode_rhs.cpp
rate_coeffs.cpp
)

FIND_PACKAGE(GSL)
IF (GSL_FOUND)
  SET (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS})
ENDIF (GSL_FOUND)

ADD_EXECUTABLE (nuclear_network ${nuclear_network_SOURCES})
TARGET_LINK_LIBRARIES(nuclear_network ${GSL_LIBRARIES})
