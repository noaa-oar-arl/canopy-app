#
# Use `FC=<compiler> ... make` to select compiler
# Use `DEBUG=1 ... make` for a debug build
# Use `NETCDF=1 ... make` to build with NetCDF using `nf-config`

# Compiler
FC ?= gfortran
ifeq ($(FC), f77)  # override possible Make default
  FC := gfortran
endif
$(info FC setting: '$(FC)')

# Default to non-debug build
DEBUG ?= 0

# Default to NetCDF build
NETCDF ?= 1

# Compile flags
$(info DEBUG setting: '$(DEBUG)')
ifeq ($(DEBUG), 1)
  FCFLAGS := -g -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
else ifeq ($(DEBUG), 0)
  FCFLAGS := -O3
else
  $(error invalid setting for DEBUG, should be 0 or 1 but is '$(DEBUG)')
endif

# NETCDF Settings here
LIBS :=
INC :=
$(info NETCDF setting: '$(NETCDF)')
ifeq ($(NETCDF), 1)
  NETCDF_FLIBS := $(shell nf-config --flibs)
  NETCDF_INC := -I$(shell nf-config --includedir)
  #
  LIBS += $(NETCDF_FLIBS)
  INC += $(NETCDF_INC)
  FCFLAGS += -DNETCDF
else ifeq ($(NETCDF), 0)
  #
else
  $(error invalid setting for NETCDF, should be 0 or 1 but is '$(NETCDF)')
endif
$(info LIBS: '$(LIBS)')
$(info INC:  '$(INC)')

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
 canopy_ncf_io_mod.o \
 canopy_check_input.o \
 canopy_read_txt.o \
 canopy_dxcalc_mod.o \
 canopy_profile_mod.o \
 canopy_wind_mod.o \
 canopy_waf_mod.o \
 canopy_phot_mod.o \
 canopy_eddy_mod.o \
 canopy_bioparm_mod.o \
 canopy_bioemi_mod.o \
 canopy_models_mod.o \
 canopy_calcs.o \
 canopy_write_txt.o \
 canopy_dealloc.o \
 canopy_app.o

ifeq ($(NETCDF), 0)
  _ncf_objs := canopy_check_input.o canopy_ncf_io_mod.o
  OBJS := $(filter-out $(_ncf_objs),$(OBJS))
endif

# Program name
PROGRAM := ../canopy

# Targets
.PHONY: all clean
all: $(PROGRAM)

$(PROGRAM): $(OBJS_APP) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $@ $(LIBS) $(INC)

%.o: %.F90
	$(FC) $(FCFLAGS) $(INC) -c $<

clean:
	rm -f *.o *.mod $(PROGRAM)
