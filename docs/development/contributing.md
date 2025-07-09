# Contributing to Canopy-App

Welcome to the Canopy-App development community! This guide will help you contribute effectively to the project.

## Getting Started

### Development Environment Setup

#### Prerequisites

```bash
# Required tools
sudo apt-get install gfortran netcdf-bin libnetcdf-dev
sudo apt-get install git make cmake
sudo apt-get install python3 python3-pip

# Optional but recommended
sudo apt-get install valgrind gdb
pip3 install pre-commit black flake8
```

#### Fork and Clone

```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR_USERNAME/canopy-app.git
cd canopy-app

# Add upstream remote
git remote add upstream https://github.com/canopy-app/canopy-app.git

# Set up pre-commit hooks
pre-commit install
```

### Development Workflow

#### 1. Create Feature Branch

```bash
# Update your fork
git checkout main
git pull upstream main
git push origin main

# Create feature branch
git checkout -b feature/your-feature-name
```

#### 2. Make Changes

```bash
# Make your changes
# Test thoroughly
# Update documentation

# Stage changes
git add .
git commit -m "feat: descriptive commit message"
```

#### 3. Submit Pull Request

```bash
# Push to your fork
git push origin feature/your-feature-name

# Create PR on GitHub
# Fill out PR template completely
```

## Code Organization

### Directory Structure

```
canopy-app/
├── src/                    # Fortran source code
│   ├── canopy_app.F90     # Main program
│   ├── *_mod.F90          # Module files
│   └── Makefile           # Build system
├── input/                 # Example input files
├── output/                # Example outputs (gitignored)
├── docs/                  # Documentation source
├── python/                # Python utilities
├── tests/                 # Test suite (future)
└── .github/               # GitHub workflows
```

### Module Organization

#### Core Modules
- `canopy_app.F90` - Main program entry point
- `canopy_init.F90` - Initialization routines
- `canopy_calcs.F90` - Main calculation driver

#### I/O Modules
- `canopy_files_mod.F90` - File management
- `canopy_ncf_io_mod.F90` - NetCDF operations
- `canopy_txt_io_mod.F90` - Text file operations

#### Physics Modules
- `canopy_canmet_mod.F90` - Meteorology
- `canopy_rad_mod.F90` - Radiation
- `canopy_bioemi_mod.F90` - Biogenic emissions
- `canopy_drydep_mod.F90` - Dry deposition

#### Support Modules
- `canopy_const_mod.F90` - Constants
- `canopy_utils_mod.F90` - Utilities
- `canopy_coord_mod.F90` - Coordinates

## Coding Standards

### Fortran Style Guide

#### Naming Conventions

```fortran
! Module names: lowercase with underscores
module canopy_physics_mod

! Subroutine/function names: lowercase with underscores
subroutine calc_wind_profile()
real function interpolate_linear()

! Variable names: descriptive, lowercase with underscores
real :: temperature_k        ! Temperature in Kelvin
real :: wind_speed_ms       ! Wind speed in m/s
integer :: num_canopy_levels

! Constants: uppercase with underscores
real, parameter :: PI = 3.14159265359
real, parameter :: GRAV_ACCEL = 9.81
```

#### Code Formatting

```fortran
! Indentation: 2 spaces (no tabs)
module example_mod
  implicit none

  ! Module variables
  real :: module_variable

  contains

  subroutine example_routine(input_var, output_var)
    ! Subroutine arguments
    real, intent(in) :: input_var
    real, intent(out) :: output_var

    ! Local variables
    integer :: i, j
    real :: temporary_value

    ! Main code block
    do i = 1, 10
      if (input_var > 0.0) then
        temporary_value = sqrt(input_var)
        output_var = temporary_value * 2.0
      else
        output_var = 0.0
      end if
    end do

  end subroutine example_routine

end module example_mod
```

#### Documentation Standards

```fortran
!> \brief Calculate wind speed profile through canopy
!!
!! This subroutine computes the vertical wind speed profile
!! through the canopy using an exponential attenuation model.
!!
!! \param[in] height_levels Vertical coordinate array (m)
!! \param[in] reference_wind Reference wind speed (m/s)
!! \param[in] reference_height Reference height (m)
!! \param[in] canopy_height Canopy top height (m)
!! \param[in] attenuation_coeff Attenuation coefficient
!! \param[out] wind_profile Computed wind speeds (m/s)
!!
!! \author Your Name
!! \date 2024-01-01
!! \version 1.0
!!
!! \see Harman and Finnigan (2007) for theoretical background
!!
subroutine calc_wind_profile(height_levels, reference_wind, &
                             reference_height, canopy_height, &
                             attenuation_coeff, wind_profile)
```

### Error Handling

#### Robust Error Checking

```fortran
subroutine safe_file_operations()
  implicit none

  integer :: iostat_val
  character(len=256) :: error_msg

  ! Check file existence
  if (.not. file_exists(input_filename)) then
    write(error_msg, '("Input file not found: ", a)') trim(input_filename)
    call abort_with_message(error_msg)
  end if

  ! Safe file opening
  open(unit=10, file=input_filename, status='old', &
       action='read', iostat=iostat_val)
  if (iostat_val /= 0) then
    write(error_msg, '("Error opening file: ", a)') trim(input_filename)
    call abort_with_message(error_msg)
  end if

  ! Always close files
  close(10)

end subroutine safe_file_operations
```

#### Input Validation

```fortran
subroutine validate_inputs()
  implicit none

  ! Check physical bounds
  if (canopy_height <= 0.0) then
    call abort_with_message("Canopy height must be positive")
  end if

  if (num_levels < 3) then
    call abort_with_message("Minimum 3 vertical levels required")
  end if

  ! Check array bounds
  if (size(temperature_array) /= num_levels) then
    call abort_with_message("Temperature array size mismatch")
  end if

  ! Warn about questionable values
  if (canopy_height > 100.0) then
    call warning_message("Canopy height > 100m is unusual")
  end if

end subroutine validate_inputs
```

<!-- ## Testing Guidelines

### Unit Testing Framework

#### Test Structure

```fortran
! test_canopy_wind.F90
program test_canopy_wind
  use canopy_wind_mod
  implicit none

  call test_exponential_profile()
  call test_above_canopy_profile()
  call test_boundary_conditions()

  write(*,*) "All wind profile tests passed!"

contains

  subroutine test_exponential_profile()
    real, parameter :: tolerance = 1.0e-6
    real :: result, expected

    ! Test known case
    call calc_wind_at_height(height=10.0, canopy_top=20.0, &
                            wind_top=5.0, attenuation=2.0, &
                            result=result)

    expected = 5.0 * exp(-2.0 * (10.0/20.0 - 1.0))

    if (abs(result - expected) > tolerance) then
      write(*,*) "FAIL: Exponential profile test"
      stop 1
    end if

  end subroutine test_exponential_profile

end program test_canopy_wind
``` -->

<!-- #### Running Tests

```bash
# Compile and run tests
cd tests
make test_canopy_wind
./test_canopy_wind

# Run all tests
make test
```

### Integration Testing

#### Full Model Tests

```bash
#!/bin/bash
# integration_test.sh

set -e  # Exit on any error

echo "Running integration tests..."

# Test 1: Basic NetCDF I/O
echo "Testing NetCDF I/O..."
cd test_data
../src/canopy_app -n test_netcdf.nml
if [ ! -f "test_output.nc" ]; then
    echo "FAIL: NetCDF output not created"
    exit 1
fi

# Test 2: Text I/O compatibility
echo "Testing text I/O..."
../src/canopy_app -n test_text.nml
if [ ! -f "test_output.txt" ]; then
    echo "FAIL: Text output not created"
    exit 1
fi

# Test 3: Regression test
echo "Testing numerical consistency..."
python3 compare_outputs.py reference_output.nc test_output.nc
if [ $? -ne 0 ]; then
    echo "FAIL: Results differ from reference"
    exit 1
fi

echo "All integration tests passed!"
```

### Performance Testing

#### Benchmarking

```fortran
! performance_timer.F90
module performance_timer
  implicit none

contains

  subroutine start_timer(timer_name)
    character(len=*), intent(in) :: timer_name
    ! Implementation using system_clock
  end subroutine start_timer

  subroutine end_timer(timer_name)
    character(len=*), intent(in) :: timer_name
    ! Implementation and reporting
  end subroutine end_timer

end module performance_timer
``` -->

## Adding New Features

### Physics Modules

#### Template for New Physics Module

```fortran
!> \file canopy_newphysics_mod.F90
!> \brief New physics process module template
!!
!! Template for implementing new physics processes in the canopy model.
!!
!! \author Your Name
!! \date 2024-01-01

!> \defgroup newphysics_mod New Physics Process Module
!! \brief Implementation of new physics process
!!
!! Detailed description of the new physics process, including:
!! - Physical background and theory
!! - Key equations and algorithms
!! - Input requirements and outputs
!! - Validation data and references
!! \{

module canopy_newphysics_mod
  use canopy_const_mod
  use canopy_utils_mod
  implicit none

  private

  ! Public interfaces
  public :: init_newphysics
  public :: calc_newphysics
  public :: finalize_newphysics

  ! Module parameters
  real, parameter :: PROCESS_CONSTANT = 1.0

  ! Module variables
  logical :: module_initialized = .false.

contains

  !> \brief Initialize new physics module
  !! \param[in] config_params Configuration parameters
  subroutine init_newphysics(config_params)
    type(config_type), intent(in) :: config_params

    ! Initialization code
    module_initialized = .true.

  end subroutine init_newphysics

  !> \brief Calculate new physics process
  !! \param[in] input_data Input meteorological data
  !! \param[out] output_data Calculated process rates
  subroutine calc_newphysics(input_data, output_data)
    type(input_type), intent(in) :: input_data
    type(output_type), intent(out) :: output_data

    ! Check initialization
    if (.not. module_initialized) then
      call abort_with_message("New physics module not initialized")
    end if

    ! Main calculations
    ! ... implementation ...

  end subroutine calc_newphysics

  !> \brief Finalize new physics module
  subroutine finalize_newphysics()

    ! Cleanup code
    module_initialized = .false.

  end subroutine finalize_newphysics

end module canopy_newphysics_mod
!! \}
```

#### Integration Steps

1. **Add module to Makefile**:
```makefile
MODULES += canopy_newphysics_mod.o
```

2. **Add namelist options**:
```fortran
! In canopy_canopts_mod.F90
logical :: opt_newphysics = .false.
```

3. **Integrate into main calculation**:
```fortran
! In canopy_calcs.F90
if (opt_newphysics) then
  call calc_newphysics(met_data, process_output)
end if
```

4. **Add output variables**:
```fortran
! In canopy_ncf_io_mod.F90
call write_var_3d(ncid, 'newphysics_rate', process_output)
```

### I/O Enhancements

#### Adding New Input Variables

```fortran
! In canopy_ncf_io_mod.F90
subroutine read_new_variable(ncid, var_name, data_array)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: var_name
  real, intent(out) :: data_array(:,:)

  integer :: varid, iostat

  ! Get variable ID
  iostat = nf90_inq_varid(ncid, var_name, varid)
  if (iostat /= nf90_noerr) then
    call handle_netcdf_error(iostat, "Variable not found: "//var_name)
  end if

  ! Read data
  iostat = nf90_get_var(ncid, varid, data_array)
  if (iostat /= nf90_noerr) then
    call handle_netcdf_error(iostat, "Error reading: "//var_name)
  end if

end subroutine read_new_variable
```

## Code Review Process

### Pull Request Guidelines

#### PR Checklist

- [ ] **Code compiles** without warnings
- [ ] **Tests pass** (unit and integration)
- [ ] **Documentation updated** (inline and user docs)
- [ ] **Performance impact** assessed
- [ ] **Backward compatibility** maintained
- [ ] **Code style** follows guidelines

#### PR Description Template

```markdown
## Description
Brief description of changes and motivation.

## Type of Change
- [ ] Bug fix (non-breaking change)
- [ ] New feature (non-breaking change)
- [ ] Breaking change (fix or feature causing existing functionality to change)
- [ ] Documentation update

## Testing
- [ ] Unit tests added/updated
- [ ] Integration tests pass
- [ ] Manual testing performed

## Performance Impact
- [ ] No performance impact
- [ ] Minor performance improvement
- [ ] Significant performance change (describe)

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review performed
- [ ] Documentation updated
- [ ] Tests added for new functionality
```

### Review Criteria

#### Technical Review

1. **Correctness**
   - Algorithm implementation
   - Mathematical accuracy
   - Boundary condition handling

2. **Performance**
   - Computational efficiency
   - Memory usage
   - Scalability considerations

3. **Maintainability**
   - Code clarity and readability
   - Modular design
   - Documentation quality

#### Science Review

1. **Physical Validity**
   - Conservation principles
   - Physical reasoning
   - Dimensional analysis

2. **Literature Basis**
   - Peer-reviewed references
   - Standard methodologies
   - Validation against observations

## Release Process

### Version Management

#### Semantic Versioning

- **MAJOR.MINOR.PATCH** (e.g., 1.2.3)
- **MAJOR**: Breaking changes
- **MINOR**: New features (backward compatible)
- **PATCH**: Bug fixes

#### Release Checklist

```bash
# 1. Update version numbers
vim src/canopy_const_mod.F90  # Update VERSION_STRING

# 2. Update CHANGELOG.md
vim CHANGELOG.md

# 3. Run full test suite
make clean && make test

# 4. Build documentation
mkdocs build

# 5. Create release tag
git tag -a v1.2.3 -m "Release version 1.2.3"
git push origin v1.2.3

# 6. Create GitHub release with binaries
```

## Getting Help

### Communication Channels

- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and ideas

### Mentoring Program

New contributors can request mentoring:
1. Create issue with "mentoring" label
2. Describe your background and interests
3. Experienced developer will be assigned
4. Regular check-ins and code reviews

## Next Steps

- **[Code Style Guide](code-style.md)** - Detailed style requirements
<!-- - **[Testing Guide](testing.md)** - Comprehensive testing procedures -->
- **[Architecture Overview](architecture.md)** - System design details
- **[API Reference](../api/overview.md)** - Code documentation
