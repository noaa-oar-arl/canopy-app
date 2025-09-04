# Code Style Guide

Comprehensive coding standards and style guidelines for Canopy-App development.

## Fortran Coding Standards

### General Principles

1. **Clarity over cleverness** - Code should be self-documenting
2. **Consistency** - Follow established patterns throughout the codebase
3. **Modularity** - Use modules to organize functionality
4. **Documentation** - All public interfaces must be documented

### File Organization

#### File Naming

```bash
# Main program
canopy_app.F90

# Module files
canopy_<functionality>_mod.F90

# Examples:
canopy_bioemi_mod.F90      # Biogenic emissions
canopy_canmet_mod.F90      # Canopy meteorology
canopy_utils_mod.F90       # Utility functions
```

#### File Structure Template

```fortran
!> \file canopy_example_mod.F90
!> \brief Brief description of module purpose
!!
!! Detailed description of the module functionality,
!! including key algorithms, assumptions, and usage.
!!
!! \author Author Name
!! \date Creation date
!! \version Version number

!> \defgroup example_mod Example Module
!! \brief Example module for demonstration
!!
!! This module demonstrates proper code organization
!! and documentation standards for the Canopy-App project.
!! \{

module canopy_example_mod

  ! Use statements (alphabetical order)
  use canopy_const_mod, only: rk => real_kind, pi, grav
  use canopy_utils_mod, only: interpolate_linear, check_bounds
  use netcdf

  ! Always include implicit none
  implicit none

  ! Module accessibility (private by default)
  private

  ! Public interfaces (explicitly list all public items)
  public :: init_example_module
  public :: calc_example_process
  public :: finalize_example_module
  public :: example_type

  ! Module parameters (uppercase with underscores)
  real(rk), parameter :: DEFAULT_VALUE = 1.0_rk
  real(rk), parameter :: MAX_ITERATIONS = 100
  character(len=*), parameter :: MODULE_NAME = 'example_mod'

  ! Module variables (lowercase with underscores)
  logical :: module_initialized = .false.
  integer :: debug_level = 0

  ! Derived types
  type :: example_type
    real(rk) :: temperature    !< Temperature (K)
    real(rk) :: pressure      !< Pressure (Pa)
    real(rk), allocatable :: profile(:)  !< Vertical profile
  end type example_type

contains

  ! Module procedures...

end module canopy_example_mod
!! \}
```

### Naming Conventions

#### Variables

```fortran
! Variables: lowercase with underscores, descriptive names
real(rk) :: air_temperature_k         ! Temperature in Kelvin
real(rk) :: wind_speed_ms            ! Wind speed in m/s
real(rk) :: leaf_area_index          ! Leaf area index (m²/m²)
integer :: num_canopy_levels         ! Number of vertical levels
logical :: use_advanced_physics      ! Physics option flag

! Arrays: plural or indicate dimensionality
real(rk) :: height_levels(nlev)      ! Height coordinates
real(rk) :: temperature_profile(:)   ! Temperature vs height
real(rk) :: emission_rates(:,:)      ! 2D emission array

! Loop indices: simple, single letters for short loops
integer :: i, j, k, n
! Descriptive names for complex loops
integer :: level_idx, species_idx, time_step
```

#### Constants

```fortran
! Constants: uppercase with underscores
real(rk), parameter :: PI = 3.14159265358979_rk
real(rk), parameter :: GRAVITATIONAL_ACCEL = 9.80665_rk
real(rk), parameter :: GAS_CONSTANT_DRY_AIR = 287.0_rk  ! J/kg/K
real(rk), parameter :: STEFAN_BOLTZMANN = 5.67e-8_rk    ! W/m²/K⁴

! Physical constants with units in comments
real(rk), parameter :: AVOGADRO_NUMBER = 6.022e23_rk    ! molecules/mol
real(rk), parameter :: BOLTZMANN_CONST = 1.381e-23_rk   ! J/K

! Model-specific constants
real(rk), parameter :: DEFAULT_CANOPY_HEIGHT = 20.0_rk  ! m
integer, parameter :: MAX_CANOPY_LEVELS = 100
integer, parameter :: MAX_SPECIES = 50
```

#### Subroutines and Functions

```fortran
! Subroutines: verb phrases, lowercase with underscores
subroutine calc_wind_profile()
subroutine read_input_data()
subroutine write_output_files()
subroutine initialize_physics_modules()

! Functions: noun phrases or questions, return type obvious
real(rk) function interpolate_linear(x, y, xi)
logical function file_exists(filename)
integer function count_valid_data_points(data_array)
real(rk) function saturation_vapor_pressure(temperature)
```

### Code Formatting

#### Indentation

```fortran
! Use 2 spaces for indentation (no tabs)
module example_mod
  implicit none

  contains

  subroutine example_routine()
    integer :: i
    real(rk) :: value

    do i = 1, 10
      if (i > 5) then
        value = real(i, rk) * 2.0_rk
      else
        value = real(i, rk)
      end if
    end do

  end subroutine example_routine

end module example_mod
```

#### Line Length

```fortran
! Maximum 80 characters per line
! Break long lines using continuation characters
subroutine long_parameter_list(temperature, pressure, humidity, &
                              wind_speed, wind_direction, &
                              solar_radiation, output_array)

! Break long expressions logically
total_emission = isoprene_emission + monoterpene_emission + &
                 sesquiterpene_emission + other_voc_emission

! Align continuation lines
if ((temperature > freezing_point) .and. &
    (pressure > minimum_pressure) .and. &
    (humidity >= 0.0_rk) .and. &
    (humidity <= 100.0_rk)) then
  ! Process valid meteorological data
end if
```

#### Spacing

```fortran
! Spaces around operators
result = a + b * c / d
logical_result = (x > 0.0_rk) .and. (y < 100.0_rk)

! Spaces after commas
call subroutine(arg1, arg2, arg3)
array = [1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk]

! No space before/after parentheses in function calls
value = sqrt(temperature * pressure)

! Space after keywords
if (condition) then
do i = 1, n
while (convergence_check)
```

### Documentation Standards

#### Module Documentation

```fortran
!> \brief Calculate biogenic emissions from vegetation
!!
!! This module implements the MEGAN (Model of Emissions of Gases
!! and Aerosols from Nature) algorithms for calculating biogenic
!! volatile organic compound emissions from vegetation.
!!
!! The module includes:
!! - Temperature and light-dependent emission factors
!! - Species-specific emission algorithms
!! - Canopy environment corrections
!! - Integration with canopy physics
!!
!! \author P.C. Campbell, NOAA ARL
!! \date 2024-01-15
!! \version 2.1
!!
!! \see Guenther et al. (2012) for algorithm details
!! \see https://bai.ess.uci.edu/megan for model documentation
!!
module canopy_bioemi_mod
```

#### Subroutine Documentation

```fortran
!> \brief Calculate isoprene emission rates for canopy layers
!!
!! Computes isoprene emission rates using the MEGAN algorithm
!! with temperature and PAR (photosynthetically active radiation)
!! dependencies. Emissions are calculated for each canopy layer
!! based on leaf area density and environmental conditions.
!!
!! \param[in] temperature Layer temperature array (K)
!! \param[in] par_radiation PAR radiation array (μmol/m²/s)
!! \param[in] leaf_area_density LAD profile (m²/m³)
!! \param[in] emission_factor Base emission factor (μg/g/h)
!! \param[out] emission_rate Calculated emission rates (μg/m³/s)
!! \param[out] status Return status (0=success, >0=error)
!!
!! \pre temperature values must be > 0 K
!! \pre par_radiation values must be >= 0
!! \pre leaf_area_density values must be >= 0
!! \post emission_rate contains valid emission rates
!!
!! \note Emission factor should be at standard conditions
!!       (30°C, 1000 μmol/m²/s PAR)
!!
!! \warning Large emission factors (>100 μg/g/h) should be verified
!!
!! \author P.C. Campbell
!! \date 2024-01-15
!! \version 1.0
!!
!! \see calc_temperature_activity() for temperature dependence
!! \see calc_light_activity() for PAR dependence
!!
subroutine calc_isoprene_emission(temperature, par_radiation, &
                                 leaf_area_density, emission_factor, &
                                 emission_rate, status)
```

#### Variable Documentation

```fortran
! In-line comments for variable declarations
real(rk) :: temperature_k          !< Air temperature (K)
real(rk) :: relative_humidity      !< Relative humidity (%)
real(rk) :: wind_speed_ms          !< Wind speed (m/s)
integer :: num_iterations          !< Number of solver iterations

! Block comments for complex variable groups
!> \defgroup met_inputs Meteorological Input Variables
!! Input meteorological variables required for canopy calculations
!! \{
real(rk) :: air_pressure_pa        !< Atmospheric pressure (Pa)
real(rk) :: solar_zenith_angle     !< Solar zenith angle (radians)
real(rk) :: cloud_fraction         !< Cloud fraction (0-1)
!! \}
```

### Error Handling

#### Error Checking Patterns

```fortran
! Input validation with descriptive messages
subroutine validate_meteorological_inputs(temp, press, humid, status)
  real(rk), intent(in) :: temp, press, humid
  integer, intent(out) :: status

  status = 0  ! Success

  ! Temperature checks
  if (temp <= 0.0_rk) then
    write(error_unit, '(a,f8.2)') 'Invalid temperature: ', temp
    status = 1
    return
  end if

  if (temp > 350.0_rk) then
    write(error_unit, '(a,f8.2)') 'Unrealistic temperature: ', temp
    status = 2
    return
  end if

  ! Pressure checks
  if (press < 50000.0_rk .or. press > 110000.0_rk) then
    write(error_unit, '(a,f10.1)') 'Invalid pressure: ', press
    status = 3
    return
  end if

  ! Humidity checks
  if (humid < 0.0_rk .or. humid > 100.0_rk) then
    write(error_unit, '(a,f8.2)') 'Invalid humidity: ', humid
    status = 4
    return
  end if

end subroutine validate_meteorological_inputs
```

#### File I/O Error Handling

```fortran
! Safe file operations with error checking
subroutine safe_file_read(filename, data_array, status)
  character(len=*), intent(in) :: filename
  real(rk), intent(out) :: data_array(:)
  integer, intent(out) :: status

  integer :: file_unit, iostat
  logical :: file_exists

  ! Check file existence
  inquire(file=filename, exist=file_exists)
  if (.not. file_exists) then
    write(error_unit, '(a)') 'File not found: ' // trim(filename)
    status = -1
    return
  end if

  ! Open file safely
  open(newunit=file_unit, file=filename, status='old', &
       action='read', iostat=iostat)
  if (iostat /= 0) then
    write(error_unit, '(a,i0)') 'Error opening file, iostat: ', iostat
    status = -2
    return
  end if

  ! Read data with error checking
  read(file_unit, *, iostat=iostat) data_array
  if (iostat /= 0) then
    write(error_unit, '(a,i0)') 'Error reading data, iostat: ', iostat
    close(file_unit)
    status = -3
    return
  end if

  close(file_unit)
  status = 0  ! Success

end subroutine safe_file_read
```

### Performance Guidelines

#### Efficient Array Operations

```fortran
! Prefer array operations over explicit loops
! Good: vectorized operation
temperature_c = temperature_k - 273.15_rk

! Less efficient: explicit loop
do i = 1, size(temperature_k)
  temperature_c(i) = temperature_k(i) - 273.15_rk
end do

! Use intrinsic functions when possible
total_emission = sum(emission_rates)
max_temperature = maxval(temperature_profile)
average_wind = sum(wind_speeds) / real(size(wind_speeds), rk)
```

#### Memory Management

```fortran
! Allocate arrays only when needed
subroutine process_canopy_data(nlevels)
  integer, intent(in) :: nlevels

  real(rk), allocatable :: work_array(:)
  integer :: alloc_stat

  ! Allocate with error checking
  allocate(work_array(nlevels), stat=alloc_stat)
  if (alloc_stat /= 0) then
    call abort_with_message('Memory allocation failed')
  end if

  ! Use array...

  ! Always deallocate
  deallocate(work_array)

end subroutine process_canopy_data
```

#### Avoiding Common Performance Issues

```fortran
! Avoid repeated function calls in loops
! Bad:
do i = 1, n
  result(i) = expensive_function(data(i))
end do

! Good: cache expensive calculations
expensive_value = expensive_function(data(1))
do i = 1, n
  if (data(i) /= data(1)) then
    expensive_value = expensive_function(data(i))
  end if
  result(i) = expensive_value
end do

! Use appropriate variable types
! Use single precision for large arrays if appropriate
real(real32) :: large_array(1000000)  ! If double precision not needed

! Avoid unnecessary type conversions
integer :: count
real(rk) :: average
average = sum(values) / real(count, rk)  ! Explicit conversion
```

## Python Coding Standards

### PEP 8 Compliance

#### Naming Conventions

```python
# Variables and functions: lowercase with underscores
air_temperature = 25.0
wind_speed_ms = 5.2

def calculate_emission_rate():
    pass

def read_meteorological_data():
    pass

# Constants: uppercase with underscores
GRAVITY = 9.81
MOLECULAR_WEIGHT_AIR = 28.97

# Classes: CapWords (PascalCase)
class CanopyModel:
    def __init__(self):
        pass

class EmissionCalculator:
    def __init__(self):
        pass
```

#### Function Documentation

```python
def calc_biogenic_emissions(temperature, par, leaf_area, emission_factor):
    """Calculate biogenic emission rates using MEGAN algorithms.

    Computes volatile organic compound emission rates from vegetation
    based on temperature, photosynthetically active radiation (PAR),
    and leaf area density.

    Args:
        temperature (np.ndarray): Air temperature array (K)
        par (np.ndarray): PAR radiation array (μmol/m²/s)
        leaf_area (np.ndarray): Leaf area density (m²/m³)
        emission_factor (float): Base emission factor (μg/g/h)

    Returns:
        np.ndarray: Emission rates (μg/m³/s)

    Raises:
        ValueError: If input arrays have different shapes
        ValueError: If temperature values are non-positive

    Example:
        >>> temp = np.array([298.15, 301.15, 295.15])
        >>> par = np.array([800.0, 1200.0, 400.0])
        >>> lad = np.array([0.5, 0.8, 0.3])
        >>> emissions = calc_biogenic_emissions(temp, par, lad, 24.0)
    """
    # Validate inputs
    if temperature.shape != par.shape != leaf_area.shape:
        raise ValueError("Input arrays must have the same shape")

    if np.any(temperature <= 0):
        raise ValueError("Temperature values must be positive")

    # Calculate emission activity factors
    gamma_t = temperature_activity_factor(temperature)
    gamma_p = par_activity_factor(par)

    # Compute emission rates
    emission_rates = emission_factor * gamma_t * gamma_p * leaf_area

    return emission_rates
```

### Code Organization

#### Module Structure

```python
"""
Canopy model utilities for data processing and analysis.

This module provides utility functions for processing meteorological
data, calculating derived variables, and interfacing with the
Fortran canopy model.

Author: Your Name
Date: 2024-01-15
Version: 1.0
"""

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from typing import Union, Tuple, Optional, Dict, Any

# Module constants
DEFAULT_CANOPY_HEIGHT = 20.0  # meters
STANDARD_TEMPERATURE = 30.0   # Celsius
STANDARD_PAR = 1000.0        # μmol/m²/s

# Type aliases for clarity
Temperature = np.ndarray
Pressure = np.ndarray
EmissionRate = np.ndarray

class CanopyDataProcessor:
    """Process and validate canopy model input/output data."""

    def __init__(self, config_file: Optional[Path] = None):
        """Initialize data processor with optional configuration."""
        self.config = self._load_config(config_file)
        self.logger = self._setup_logging()

    def process_meteorological_data(self,
                                  input_file: Path,
                                  output_file: Path) -> Dict[str, Any]:
        """Process meteorological input data for canopy model."""
        # Implementation...
        pass
```

## Build System Standards

### Makefile Organization

```makefile
# Canopy-App Makefile
# Author: P.C. Campbell
# Date: 2024-01-15

# Compiler settings
FC = gfortran
FFLAGS = -O2 -g -Wall -Wextra -fcheck=bounds
LDFLAGS =
LIBS = -lnetcdff -lnetcdf

# Directories
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# Source files (in dependency order)
SOURCES = canopy_const_mod.F90 \
          canopy_utils_mod.F90 \
          canopy_coord_mod.F90 \
          canopy_files_mod.F90 \
          canopy_bioemi_mod.F90 \
          canopy_canmet_mod.F90 \
          canopy_rad_mod.F90 \
          canopy_drydep_mod.F90 \
          canopy_phot_mod.F90 \
          canopy_ncf_io_mod.F90 \
          canopy_txt_io_mod.F90 \
          canopy_init.F90 \
          canopy_alloc.F90 \
          canopy_dealloc.F90 \
          canopy_calcs.F90 \
          canopy_app.F90

# Object files
OBJECTS = $(SOURCES:%.F90=$(OBJDIR)/%.o)

# Target executable
TARGET = $(BINDIR)/canopy_app

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(OBJECTS) | $(BINDIR)
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

# Compile source files
$(OBJDIR)/%.o: $(SRCDIR)/%.F90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Create directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# Clean build files
clean:
	rm -rf $(OBJDIR) $(BINDIR)

# Development targets
debug: FFLAGS += -g -O0 -fbacktrace -ffpe-trap=invalid,zero,overflow
debug: $(TARGET)

profile: FFLAGS += -pg
profile: $(TARGET)

# Testing targets
test: $(TARGET)
	cd tests && ./run_tests.sh

# Installation
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

.PHONY: all clean debug profile test install
```

## Documentation Standards

### README Structure

```markdown
# Canopy-App

Brief description of the project and its purpose.

## Quick Start

```bash
# Clone and build
git clone https://github.com/canopy-app/canopy-app.git
cd canopy-app
make

# Run example
cd input
../bin/canopy_app
```

## Installation

Detailed installation instructions...

## Usage

Basic usage examples...

## Documentation

Link to full documentation...

## Contributing

How to contribute to the project...

## License

License information...
```

### Changelog Format

```markdown
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- New biogenic emission species
- Enhanced radiation calculations

### Changed
- Improved numerical stability in wind profile calculations
- Updated documentation format

### Fixed
- Memory leak in NetCDF I/O routines
- Boundary condition handling in dry deposition

## [1.2.0] - 2024-01-15

### Added
- Complete Doxygen documentation
- MkDocs documentation site
- Python analysis utilities

### Changed
- Refactored emission calculation algorithms
- Improved error handling throughout codebase

### Fixed
- Division by zero in stability calculations
- Incorrect units in photolysis rate output
```

## Best Practices Summary

### Code Quality Checklist

- [ ] **Meaningful names** for all variables and functions
- [ ] **Consistent indentation** (2 spaces)
- [ ] **Comprehensive documentation** for all public interfaces
- [ ] **Error checking** for all I/O operations
- [ ] **Input validation** for all subroutines
- [ ] **Memory management** (proper allocation/deallocation)
- [ ] **Performance considerations** (avoid unnecessary operations)
- [ ] **Modular design** (clear separation of concerns)

### Review Guidelines

1. **Readability**: Can another developer understand the code?
2. **Correctness**: Does the code do what it's supposed to do?
3. **Performance**: Are there any obvious inefficiencies?
4. **Maintainability**: Will this code be easy to modify later?
5. **Documentation**: Is the code adequately documented?

## Next Steps

- **[Testing Guide](testing.md)** - Comprehensive testing procedures
- **[Contributing Guide](contributing.md)** - Development workflow
- **[Architecture Overview](architecture.md)** - System design
- **[API Reference](../api/overview.md)** - Code documentation
