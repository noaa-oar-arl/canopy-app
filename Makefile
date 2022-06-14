# compiler
FC := gfortran

# compile flags
FCFLAGS = -g -c -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
# link flags
FLFLAGS =

# source files and objects
SRCS = canopy_wind_mod.F90 canopy_utils_mod.F90 canopy_driver.F90

# program name
PROGRAM = canopy

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FLFLAGS) -o $@ $^

%.mod: %.F90
	$(FC) $(FCFLAGS) -o $@ $<

%.o: %.F90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod $(PROGRAM)
