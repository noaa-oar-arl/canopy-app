#
# Use `DEBUG=1 make` for a debug build

# Compiler
FC ?= gfortran
ifeq ($(FC), f77)  # override possible Make default
  FC := gfortran
endif
$(info FC setting: '$(FC)')

# Default to non-debug build
DEBUG ?= 0

#Default to non-NetCDF build (required for now)
#NETCDF ?= 0

# Compile flags
ifeq ($(DEBUG), 1)
  FCFLAGS := -g -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
else ifeq ($(DEBUG), 0)
  FCFLAGS := -O3
else
  $(error invalid setting for DEBUG, should be 0 or 1 but is '$(DEBUG)')
endif
$(info DEBUG setting: '$(DEBUG)')

#ifeq ($(NETCDF), 1) to use NetCDF as option (required for building now, and set below)
  #Put NETCDF Settings here
  NETCDFF :=  /opt/sw/spack/apps/linux-centos8-cascadelake/gcc-9.3.0-openmpi-4.0.4/netcdf-c-4.7.4-vh
  LIBS   := -L$(NETCDFF)/lib -lnetcdf -lnetcdff -fopenmp
  FFLAGS := -I$(NETCDFF)/include
#  OBJS_APP := \
#	  canopy_ncf_io_mod.o \
#	  canopy_read_ncf.o
#else ifeq ($(NETCDF), 0)
#  FFLAGS :=
#  OBJS_APP :=
#else
#  $(error invalid setting for NETCDF, should be 0 or 1 but is '$(NETCDF)')
#endif
#$(info NETCDF setting: '$(NETCDF)')

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
 canopy_check_input.o \
 canopy_txt_io_mod.o \
 canopy_ncf_io_mod.o \
 canopy_read_txt.o \
 canopy_read_ncf.o \
 canopy_dxcalc_mod.o \
 canopy_profile_mod.o \
 canopy_wind_mod.o \
 canopy_waf_mod.o \
 canopy_phot_mod.o \
 canopy_eddy_mod.o \
 canopy_calcs.o \
 canopy_outncf_alloc.o \
 canopy_outncf_init.o \
 canopy_outncfglobal.o \
 canopy_write_txt.o \
 canopy_write_ncf.o \
 canopy_dealloc.o \
 canopy_app.o

# Program name
PROGRAM := canopy

# Targets
.PHONY: all clean
all: $(PROGRAM)

$(PROGRAM): $(OBJS_APP) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $@ $(FFLAGS) $(LIBS)

%.o: %.F90
	$(FC) $(FCFLAGS) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(PROGRAM)
