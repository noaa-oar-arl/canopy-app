#
# Use `DEBUG=1 make` for a debug build
#

# Compiler
FC ?= gfortran
ifeq ($(FC), f77)  # override possible Make default
  FC := gfortran
endif
$(info FC setting: '$(FC)')

# Default to non-debug build
DEBUG ?= 0

# Compile flags
ifeq ($(DEBUG), 1)
  FCFLAGS := -g -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
else ifeq ($(DEBUG), 0)
  FCFLAGS := -O3
else
  $(error invalid setting for DEBUG, should be 0 or 1 but is '$(DEBUG)')
endif
$(info DEBUG setting: '$(DEBUG)')

# Link flags
FLFLAGS :=

# Source objects
OBJS := canopy_files_mod.o canopy_const_mod.o canopy_parm_mod.o canopy_utils_mod.o canopy_waf_mod.o canopy_wind_mod.o canopy_readnml.o canopy_driver.o

# Program name
PROGRAM := canopy

# Targets
.PHONY: all clean
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $@ $(FLFLAGS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod $(PROGRAM)
