CXX = g++
CFLAGS = -c -O3
GSL_DIR = /home/friesen
GSL_INC = -I $(GSL_DIR)/include
GSL_LIB = -L $(GSL_DIR)/lib -lgsl -lgslcblas
LD = g++
LDFLAGS = -O3
SOURCES = \
	rate_coeffs.cpp \
	ode_rhs.cpp \
	jacobian.cpp \
	main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = run

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(LD) $(LDFLAGS) $(GSL_LIB) $(OBJECTS) -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(GSL_INC) $< -o $@

clean:
	rm -rf $(EXECUTABLE) *.o
