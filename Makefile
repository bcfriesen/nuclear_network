CC = gcc
CFLAGS = -c -O3
LD = gcc
LDFLAGS = -O3
SOURCES = \
          jacobian.cpp \
	  ode_rhs.cpp \
	  rate_coeffs.cpp \
	  main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = run

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.c
	$(CC) $(CCFLAGS) $< -o $@

clean:
	rm -rf $(EXECUTABLE) *.o
