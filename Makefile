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
OBJS :=\
 canopy_const_mod.o \
 canopy_coord_mod.o \
 canopy_canopts_mod.o \
 canopy_canmet_mod.o \
 canopy_canvars_mod.o \
 canopy_utils_mod.o \
 canopy_files_mod.o \
 canopy_readnml.o \
 canopy_alloc.o \
 canopy_init.o \
 canopy_txt_io_mod.o \
 canopy_read_txt.o \
 canopy_driver_mod.o \
 canopy_wind_mod.o \
 canopy_waf_mod.o \
 canopy_phot_mod.o \
 canopy_eddy_mod.o \
 canopy_calcs.o \
 canopy_write_txt.o \
 canopy_dealloc.o \
 canopy_app.o

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
