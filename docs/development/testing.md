# Testing

This document describes the testing framework and procedures for the Canopy-App model, including unit tests, integration tests, and validation approaches.

## Testing Framework

### Overview

The Canopy-App uses a comprehensive testing strategy:

1. **Unit Tests**: Test individual modules and subroutines
2. **Integration Tests**: Test component interactions
3. **System Tests**: Test complete model runs
4. **Validation Tests**: Compare against observations
5. **Performance Tests**: Monitor computational performance
6. **Regression Tests**: Ensure changes don't break functionality

### Test Organization

```
test/
├── unit/                    # Unit tests for individual modules
│   ├── test_canopy_rad.F90
│   ├── test_canopy_phot.F90
│   ├── test_canopy_wind.F90
│   └── ...
├── integration/             # Integration tests
│   ├── test_full_physics.F90
│   ├── test_io_chain.F90
│   └── ...
├── system/                  # Full system tests
│   ├── basic_run/
│   ├── chemistry_run/
│   └── ...
├── validation/              # Validation against observations
│   ├── fluxnet_sites/
│   ├── field_campaigns/
│   └── ...
├── performance/             # Performance benchmarks
│   ├── scaling_tests/
│   └── memory_tests/
└── framework/               # Testing utilities
    ├── test_framework.F90
    └── test_utilities.F90
```

## Unit Testing

### Test Framework

#### Test Module Template

```fortran
!> \file test_bioemi_mod.F90
!! \brief Unit tests for biogenic emission module
!! \author Test Author
!! \date 2024-07-15

program test_bioemi_mod
    use canopy_bioemi_mod
    use test_utils_mod
    use assert_mod
    implicit none

    ! Test statistics
    integer :: total_tests = 0
    integer :: passed_tests = 0
    integer :: failed_tests = 0

    call print_test_header('Biogenic Emission Module Tests')

    ! Run test suites
    call test_isoprene_calculations()
    call test_monoterpene_calculations()
    call test_temperature_dependencies()
    call test_light_dependencies()
    call test_boundary_conditions()
    call test_error_handling()

    call print_test_summary(total_tests, passed_tests, failed_tests)

    ! Exit with error code if any tests failed
    if (failed_tests > 0) stop 1

contains

    subroutine test_isoprene_calculations()
        call print_test_group('Isoprene Emission Calculations')

        call test_standard_conditions()
        call test_high_temperature()
        call test_low_light()
        call test_zero_lai()
    end subroutine

    subroutine test_standard_conditions()
        real(kind=real64) :: emission_rate, expected, tolerance
        integer :: status

        total_tests = total_tests + 1

        ! Standard conditions: 30°C, 1000 μmol/m²/s PAR, LAI=4
        call calc_isoprene_emission(temperature=303.15_real64, &
                                   par=1000.0_real64, &
                                   lai=4.0_real64, &
                                   emission_rate=emission_rate, &
                                   status=status)

        ! Expected value from Guenther et al. (2012)
        expected = 24.5_real64  ! μg/m²/s
        tolerance = 0.5_real64

        call assert_equal_real(emission_rate, expected, tolerance, &
                              'Standard conditions isoprene emission')
        call assert_equal_int(status, 0, 'Status code should be success')

        if (abs(emission_rate - expected) <= tolerance .and. status == 0) then
            passed_tests = passed_tests + 1
            call print_test_result('PASS', 'Standard conditions')
        else
            failed_tests = failed_tests + 1
            call print_test_result('FAIL', 'Standard conditions')
            print *, '  Expected:', expected, ' Got:', emission_rate
        end if
    end subroutine

    subroutine test_temperature_dependencies()
        call print_test_group('Temperature Dependencies')

        ! Test Q10 behavior
        call test_q10_temperature_response()
        call test_extreme_temperatures()
        call test_temperature_interpolation()
    end subroutine

    subroutine test_q10_temperature_response()
        real(kind=real64) :: emission_low, emission_high, q10_calculated
        real(kind=real64) :: temp_low, temp_high, expected_q10
        integer :: status

        total_tests = total_tests + 1

        temp_low = 293.15_real64   ! 20°C
        temp_high = 303.15_real64  ! 30°C

        ! Calculate emissions at two temperatures
        call calc_isoprene_emission(temperature=temp_low, &
                                   par=1000.0_real64, &
                                   lai=4.0_real64, &
                                   emission_rate=emission_low, &
                                   status=status)

        call calc_isoprene_emission(temperature=temp_high, &
                                   par=1000.0_real64, &
                                   lai=4.0_real64, &
                                   emission_rate=emission_high, &
                                   status=status)

        ! Calculate Q10 value
        q10_calculated = (emission_high / emission_low) ** (10.0 / (temp_high - temp_low))
        expected_q10 = 2.3_real64  ! Expected Q10 for isoprene

        call assert_equal_real(q10_calculated, expected_q10, 0.2_real64, &
                              'Q10 temperature response')

        if (abs(q10_calculated - expected_q10) <= 0.2_real64) then
            passed_tests = passed_tests + 1
            call print_test_result('PASS', 'Q10 temperature response')
        else
            failed_tests = failed_tests + 1
            call print_test_result('FAIL', 'Q10 temperature response')
            print *, '  Expected Q10:', expected_q10, ' Got:', q10_calculated
        end if
    end subroutine

    subroutine test_boundary_conditions()
        call print_test_group('Boundary Conditions')

        call test_zero_inputs()
        call test_negative_inputs()
        call test_extreme_values()
    end subroutine

    subroutine test_zero_inputs()
        real(kind=real64) :: emission_rate
        integer :: status

        total_tests = total_tests + 1

        ! Test with zero LAI
        call calc_isoprene_emission(temperature=303.15_real64, &
                                   par=1000.0_real64, &
                                   lai=0.0_real64, &
                                   emission_rate=emission_rate, &
                                   status=status)

        call assert_equal_real(emission_rate, 0.0_real64, 1.0e-10_real64, &
                              'Zero LAI should give zero emission')

        if (abs(emission_rate) <= 1.0e-10_real64) then
            passed_tests = passed_tests + 1
            call print_test_result('PASS', 'Zero LAI test')
        else
            failed_tests = failed_tests + 1
            call print_test_result('FAIL', 'Zero LAI test')
        end if
    end subroutine

    subroutine test_error_handling()
        call print_test_group('Error Handling')

        call test_invalid_temperature()
        call test_negative_par()
        call test_invalid_lai()
    end subroutine

    subroutine test_invalid_temperature()
        real(kind=real64) :: emission_rate
        integer :: status

        total_tests = total_tests + 1

        ! Test with invalid (negative) temperature
        call calc_isoprene_emission(temperature=-10.0_real64, &
                                   par=1000.0_real64, &
                                   lai=4.0_real64, &
                                   emission_rate=emission_rate, &
                                   status=status)

        call assert_not_equal_int(status, 0, 'Invalid temperature should return error')

        if (status /= 0) then
            passed_tests = passed_tests + 1
            call print_test_result('PASS', 'Invalid temperature error handling')
        else
            failed_tests = failed_tests + 1
            call print_test_result('FAIL', 'Invalid temperature error handling')
        end if
    end subroutine

end program test_bioemi_mod
```

### Assertion Module

```fortran
!> \file assert_mod.F90
!! \brief Assertion utilities for unit testing

module assert_mod
    implicit none

    private
    public :: assert_equal_real, assert_equal_int, assert_not_equal_int
    public :: assert_true, assert_false, assert_array_equal

    integer, parameter :: real64 = selected_real_kind(15, 307)

contains

    subroutine assert_equal_real(actual, expected, tolerance, message)
        real(kind=real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, 'ASSERTION FAILED: ', trim(message)
            print *, '  Expected: ', expected
            print *, '  Actual:   ', actual
            print *, '  Tolerance:', tolerance
            print *, '  Difference:', abs(actual - expected)
        end if
    end subroutine

    subroutine assert_equal_int(actual, expected, message)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: message

        if (actual /= expected) then
            print *, 'ASSERTION FAILED: ', trim(message)
            print *, '  Expected: ', expected
            print *, '  Actual:   ', actual
        end if
    end subroutine

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            print *, 'ASSERTION FAILED: ', trim(message)
            print *, '  Expected: TRUE'
            print *, '  Actual:   FALSE'
        end if
    end subroutine

    subroutine assert_array_equal(actual, expected, tolerance, message)
        real(kind=real64), intent(in) :: actual(:), expected(:)
        real(kind=real64), intent(in) :: tolerance
        character(len=*), intent(in) :: message

        integer :: i, n
        logical :: arrays_equal

        n = size(actual)
        if (size(expected) /= n) then
            print *, 'ASSERTION FAILED: ', trim(message)
            print *, '  Array sizes differ: ', n, ' vs ', size(expected)
            return
        end if

        arrays_equal = .true.
        do i = 1, n
            if (abs(actual(i) - expected(i)) > tolerance) then
                arrays_equal = .false.
                exit
            end if
        end do

        if (.not. arrays_equal) then
            print *, 'ASSERTION FAILED: ', trim(message)
            print *, '  Arrays differ at element ', i
            print *, '  Expected: ', expected(i)
            print *, '  Actual:   ', actual(i)
        end if
    end subroutine

end module assert_mod
```

### Test Utilities Module

```fortran
!> \file test_utils_mod.F90
!! \brief Common utilities for testing

module test_utils_mod
    implicit none

    private
    public :: print_test_header, print_test_summary, print_test_group
    public :: print_test_result, load_test_data, compare_files

contains

    subroutine print_test_header(test_name)
        character(len=*), intent(in) :: test_name

        print *, repeat('=', 60)
        print *, trim(test_name)
        print *, repeat('=', 60)
    end subroutine

    subroutine print_test_summary(total, passed, failed)
        integer, intent(in) :: total, passed, failed

        print *, repeat('-', 60)
        print *, 'TEST SUMMARY:'
        print *, '  Total tests:  ', total
        print *, '  Passed tests: ', passed
        print *, '  Failed tests: ', failed
        if (failed == 0) then
            print *, '  Result: ALL TESTS PASSED'
        else
            print *, '  Result: SOME TESTS FAILED'
        end if
        print *, repeat('-', 60)
    end subroutine

    subroutine print_test_group(group_name)
        character(len=*), intent(in) :: group_name

        print *, ''
        print *, trim(group_name), ':'
        print *, repeat('-', len_trim(group_name) + 1)
    end subroutine

    subroutine print_test_result(result, test_name)
        character(len=*), intent(in) :: result, test_name

        print '(a,1x,a)', '  ' // result // ':', trim(test_name)
    end subroutine

end module test_utils_mod
```

## Integration Testing

### Full Simulation Tests

```fortran
!> \file test_full_simulation.F90
!! \brief Integration test for complete model simulation

program test_full_simulation
    use canopy_app_mod
    use test_utils_mod
    implicit none

    call print_test_header('Full Simulation Integration Tests')

    call test_basic_simulation()
    call test_netcdf_io()
    call test_text_io()
    call test_multi_day_simulation()

contains

    subroutine test_basic_simulation()
        character(len=256) :: input_file, output_file, reference_file
        integer :: status
        logical :: files_match

        call print_test_group('Basic Simulation Test')

        input_file = '../data/input/test_basic_input.nc'
        output_file = 'test_basic_output.nc'
        reference_file = '../data/reference/basic_output_reference.nc'

        ! Run simulation
        call run_canopy_simulation(input_file, output_file, status)

        if (status /= 0) then
            print *, 'FAIL: Simulation returned error code:', status
            return
        end if

        ! Compare output with reference
        call compare_netcdf_files(output_file, reference_file, files_match)

        if (files_match) then
            print *, 'PASS: Basic simulation test'
        else
            print *, 'FAIL: Output differs from reference'
        end if
    end subroutine

end program test_full_simulation
```

### I/O Testing

```bash
#!/bin/bash
# test_io_operations.sh

echo "Testing I/O operations..."

# Test NetCDF input/output
echo "1. Testing NetCDF I/O..."
./test_netcdf_io
if [ $? -ne 0 ]; then
    echo "NetCDF I/O test failed"
    exit 1
fi

# Test text file I/O
echo "2. Testing text file I/O..."
./test_text_io
if [ $? -ne 0 ]; then
    echo "Text I/O test failed"
    exit 1
fi

# Test error handling
echo "3. Testing error handling..."
./test_io_errors
if [ $? -ne 0 ]; then
    echo "I/O error handling test failed"
    exit 1
fi

echo "All I/O tests passed!"
```

## Performance Testing

### Benchmark Framework

```fortran
!> \file benchmark_bioemi.F90
!! \brief Performance benchmark for biogenic emission module

program benchmark_bioemi
    use canopy_bioemi_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    integer, parameter :: nlevs = 50, ntimes = 1000
    real(kind=real64) :: temperature(nlevs), par(nlevs), lai(nlevs)
    real(kind=real64) :: emission_rates(nlevs)
    real(kind=real64) :: start_time, end_time, elapsed_time
    integer :: i, status

    ! Initialize test data
    call initialize_test_data()

    ! Benchmark emission calculations
    call cpu_time(start_time)

    do i = 1, ntimes
        call calc_biogenic_emissions(temperature, par, lai, &
                                    emission_rates, status)
    end do

    call cpu_time(end_time)
    elapsed_time = end_time - start_time

    ! Report results
    print *, 'Biogenic Emission Benchmark Results:'
    print *, '  Number of levels:', nlevs
    print *, '  Number of iterations:', ntimes
    print *, '  Total time (s):', elapsed_time
    print *, '  Time per call (ms):', (elapsed_time / ntimes) * 1000.0
    print *, '  Calls per second:', ntimes / elapsed_time

contains

    subroutine initialize_test_data()
        integer :: k

        do k = 1, nlevs
            temperature(k) = 298.15 + 5.0 * sin(real(k) / real(nlevs) * 3.14159)
            par(k) = 800.0 + 200.0 * cos(real(k) / real(nlevs) * 3.14159)
            lai(k) = 0.5 + 0.3 * exp(-real(k) / 10.0)
        end do
    end subroutine

end program benchmark_bioemi
```

### Performance Regression Testing

```bash
#!/bin/bash
# performance_regression.sh

echo "Running performance regression tests..."

# Create directories for results
mkdir -p performance_results/current
mkdir -p performance_results/baseline

# Run current benchmarks
echo "Running current benchmarks..."
./benchmark_bioemi > performance_results/current/bioemi.txt
./benchmark_radiation > performance_results/current/radiation.txt
./benchmark_canmet > performance_results/current/canmet.txt

# Compare with baseline (if exists)
if [ -d "performance_results/baseline" ]; then
    echo "Comparing with baseline performance..."
    python3 scripts/compare_performance.py \
        performance_results/baseline/ \
        performance_results/current/

    if [ $? -ne 0 ]; then
        echo "Performance regression detected!"
        exit 1
    fi
else
    echo "No baseline found. Current results saved as new baseline."
    cp -r performance_results/current/* performance_results/baseline/
fi

echo "Performance tests completed successfully."
```

## Scientific Validation Testing

### Field Data Comparison

```python
#!/usr/bin/env python3
"""
Validate model output against field observations.
"""

import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def validate_against_flux_tower(model_file, obs_file):
    """Compare model output with flux tower observations."""

    # Read model output
    model_data = nc.Dataset(model_file, 'r')
    model_time = model_data.variables['time'][:]
    model_isoprene = model_data.variables['isoprene_flux'][:]
    model_sensible_heat = model_data.variables['sensible_heat_flux'][:]

    # Read observations
    obs_data = pd.read_csv(obs_file, parse_dates=['datetime'])
    obs_time = obs_data['datetime']
    obs_isoprene = obs_data['isoprene_flux']
    obs_sensible_heat = obs_data['sensible_heat_flux']

    # Align time series (simplified - needs proper implementation)
    # ... time alignment code ...

    # Calculate statistics
    iso_corr, iso_p = stats.pearsonr(model_isoprene, obs_isoprene)
    iso_rmse = np.sqrt(np.mean((model_isoprene - obs_isoprene)**2))
    iso_bias = np.mean(model_isoprene - obs_isoprene)

    heat_corr, heat_p = stats.pearsonr(model_sensible_heat, obs_sensible_heat)
    heat_rmse = np.sqrt(np.mean((model_sensible_heat - obs_sensible_heat)**2))
    heat_bias = np.mean(model_sensible_heat - obs_sensible_heat)

    # Print validation results
    print("Validation Results:")
    print(f"Isoprene Flux:")
    print(f"  Correlation: {iso_corr:.3f} (p={iso_p:.3f})")
    print(f"  RMSE: {iso_rmse:.3f}")
    print(f"  Bias: {iso_bias:.3f}")

    print(f"Sensible Heat Flux:")
    print(f"  Correlation: {heat_corr:.3f} (p={heat_p:.3f})")
    print(f"  RMSE: {heat_rmse:.3f}")
    print(f"  Bias: {heat_bias:.3f}")

    # Create validation plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Isoprene scatter plot
    axes[0].scatter(obs_isoprene, model_isoprene, alpha=0.6)
    axes[0].plot([obs_isoprene.min(), obs_isoprene.max()],
                 [obs_isoprene.min(), obs_isoprene.max()], 'r--', lw=2)
    axes[0].set_xlabel('Observed Isoprene Flux')
    axes[0].set_ylabel('Modeled Isoprene Flux')
    axes[0].set_title(f'Isoprene (R={iso_corr:.3f})')

    # Sensible heat scatter plot
    axes[1].scatter(obs_sensible_heat, model_sensible_heat, alpha=0.6)
    axes[1].plot([obs_sensible_heat.min(), obs_sensible_heat.max()],
                 [obs_sensible_heat.min(), obs_sensible_heat.max()], 'r--', lw=2)
    axes[1].set_xlabel('Observed Sensible Heat Flux')
    axes[1].set_ylabel('Modeled Sensible Heat Flux')
    axes[1].set_title(f'Sensible Heat (R={heat_corr:.3f})')

    plt.tight_layout()
    plt.savefig('validation_results.png')

    # Return validation metrics
    return {
        'isoprene_correlation': iso_corr,
        'isoprene_rmse': iso_rmse,
        'heat_correlation': heat_corr,
        'heat_rmse': heat_rmse
    }

if __name__ == '__main__':
    import sys

    if len(sys.argv) != 3:
        print("Usage: python3 validate_field_data.py <model_file> <obs_file>")
        sys.exit(1)

    model_file = sys.argv[1]
    obs_file = sys.argv[2]

    metrics = validate_against_flux_tower(model_file, obs_file)

    # Check validation criteria
    if (metrics['isoprene_correlation'] < 0.7 or
        metrics['heat_correlation'] < 0.8):
        print("Validation criteria not met!")
        sys.exit(1)

    print("Validation successful!")
```

### Chamber Study Validation

```fortran
!> \file test_chamber_validation.F90
!! \brief Validate against environmental chamber studies

program test_chamber_validation
    use canopy_bioemi_mod
    implicit none

    ! Validate against Carter et al. (2000) chamber study
    call validate_carter_2000()

    ! Validate against Guenther et al. (1991) measurements
    call validate_guenther_1991()

contains

    subroutine validate_carter_2000()
        ! Implementation of chamber study validation
        real(kind=real64) :: temperatures(5) = [298.0, 303.0, 308.0, 313.0, 318.0]
        real(kind=real64) :: observed_emissions(5) = [12.5, 24.8, 45.2, 78.1, 125.6]
        real(kind=real64) :: modeled_emission, tolerance
        integer :: i, status, passed_tests

        passed_tests = 0
        tolerance = 10.0  ! 10% tolerance

        do i = 1, 5
            call calc_isoprene_emission(temperature=temperatures(i), &
                                       par=1000.0_real64, &
                                       lai=1.0_real64, &
                                       emission_rate=modeled_emission, &
                                       status=status)

            if (abs(modeled_emission - observed_emissions(i)) / observed_emissions(i) < 0.1) then
                passed_tests = passed_tests + 1
            end if
        end do

        print *, 'Carter et al. (2000) validation:'
        print *, '  Tests passed:', passed_tests, '/', 5

        if (passed_tests >= 4) then
            print *, '  Result: PASS'
        else
            print *, '  Result: FAIL'
        end if
    end subroutine

end program test_chamber_validation
```

## Continuous Integration

### GitHub Actions Workflow

```yaml
# .github/workflows/tests.yml
name: Comprehensive Testing

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main ]

jobs:
  unit-tests:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        compiler: [gfortran-9, gfortran-10]

    steps:
    - uses: actions/checkout@v3

    - name: Install dependencies
      run: |
        sudo apt-get update
        sudo apt-get install gfortran netcdf-bin libnetcdf-dev
        sudo apt-get install python3-pip
        pip3 install numpy scipy matplotlib netCDF4

    - name: Compile tests
      run: |
        cd tests/unit
        make FC=${{ matrix.compiler }}

    - name: Run unit tests
      run: |
        cd tests/unit
        ./run_unit_tests.sh

    - name: Upload test results
      uses: actions/upload-artifact@v3
      if: always()
      with:
        name: unit-test-results-${{ matrix.compiler }}
        path: tests/unit/results/

  integration-tests:
    needs: unit-tests
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v3

    - name: Setup environment
      run: |
        sudo apt-get update
        sudo apt-get install gfortran netcdf-bin libnetcdf-dev

    - name: Compile application
      run: |
        cd src
        make

    - name: Run integration tests
      run: |
        cd tests/integration
        ./run_integration_tests.sh

    - name: Validate against field data
      run: |
        cd tests/validation
        python3 validate_field_data.py ../data/model_output.nc ../data/observations.csv

  performance-tests:
    needs: unit-tests
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v3

    - name: Setup environment
      run: |
        sudo apt-get update
        sudo apt-get install gfortran netcdf-bin libnetcdf-dev

    - name: Compile benchmarks
      run: |
        cd tests/performance
        make

    - name: Run performance tests
      run: |
        cd tests/performance
        ./run_benchmarks.sh

    - name: Check for performance regression
      run: |
        cd tests/performance
        python3 check_regression.py baseline/ current/
```

### Test Automation Scripts

```bash
#!/bin/bash
# run_all_tests.sh - Master test script

set -e  # Exit on any error

echo "Running Canopy-App test suite..."

# Set test environment
export CANOPY_TEST_DATA_DIR="$(pwd)/tests/data"
export OMP_NUM_THREADS=1  # Ensure reproducible results

# Unit tests
echo "1. Running unit tests..."
cd tests/unit
make clean && make
./run_unit_tests.sh
cd ../..

# Integration tests
echo "2. Running integration tests..."
cd tests/integration
./run_integration_tests.sh
cd ../..

# Performance tests
echo "3. Running performance tests..."
cd tests/performance
./run_benchmarks.sh
cd ../..

# Scientific validation
echo "4. Running scientific validation..."
cd tests/validation
python3 run_validation_suite.py
cd ../..

# Memory leak checks (if valgrind available)
if command -v valgrind &> /dev/null; then
    echo "5. Running memory leak checks..."
    cd tests/memory
    ./run_memory_tests.sh
    cd ../..
fi

echo "All tests completed successfully!"
```

## Test Data Management

### Test Data Organization

```
tests/data/
├── input/
│   ├── simple_case.nc         # Basic test case
│   ├── complex_case.nc        # Multi-physics test case
│   ├── boundary_test.nc       # Edge cases
│   └── error_test.nc          # Invalid data for error testing
├── reference/
│   ├── simple_output.nc       # Expected output for simple case
│   ├── complex_output.nc      # Expected output for complex case
│   └── benchmark_results.txt  # Performance benchmarks
└── validation/
    ├── flux_tower_data.csv    # Field observations
    ├── chamber_studies.csv    # Laboratory measurements
    └── model_intercomparison/ # Other model results
```

### Data Generation Scripts

```python
#!/usr/bin/env python3
"""
Generate synthetic test data for unit and integration tests.
"""

import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta

def create_simple_test_case():
    """Create a simple test case with known analytical solutions."""

    # Define dimensions
    nlevs = 20
    ntimes = 24

    # Create NetCDF file
    with nc.Dataset('tests/data/input/simple_case.nc', 'w') as dataset:
        # Dimensions
        dataset.createDimension('levels', nlevs)
        dataset.createDimension('time', ntimes)

        # Coordinate variables
        levels = dataset.createVariable('levels', 'f8', ('levels',))
        time = dataset.createVariable('time', 'f8', ('time',))

        # Meteorological variables
        temperature = dataset.createVariable('temperature', 'f8', ('time', 'levels'))
        wind_speed = dataset.createVariable('wind_speed', 'f8', ('time', 'levels'))
        par = dataset.createVariable('par', 'f8', ('time', 'levels'))

        # Fill coordinate variables
        levels[:] = np.linspace(0, 25, nlevs)  # 0 to 25 meters
        time[:] = np.arange(ntimes)  # Hours

        # Fill meteorological data with simple patterns
        for t in range(ntimes):
            # Diurnal temperature cycle
            temp_surface = 288 + 8 * np.sin(2 * np.pi * t / 24)
            temperature[t, :] = temp_surface + 0.1 * levels[:]

            # Wind profile (log profile above canopy, exponential in canopy)
            wind_ref = 3.0 + 2.0 * np.sin(2 * np.pi * t / 24)
            for k in range(nlevs):
                if levels[k] > 20:  # Above canopy
                    wind_speed[t, k] = wind_ref * np.log((levels[k] - 14) / 2) / np.log(6 / 2)
                else:  # In canopy
                    wind_speed[t, k] = wind_ref * np.exp(-2 * (20 - levels[k]) / 20)

            # PAR with simple diurnal cycle
            if 6 <= t <= 18:  # Daytime
                par_top = 1500 * np.sin(np.pi * (t - 6) / 12)
                for k in range(nlevs):
                    # Exponential attenuation in canopy
                    lai_above = 4.0 * np.exp(-(25 - levels[k]) / 10)
                    par[t, k] = par_top * np.exp(-0.5 * lai_above)
            else:  # Nighttime
                par[t, :] = 0.0

        # Add attributes
        dataset.title = "Simple test case for Canopy-App"
        dataset.created = datetime.now().isoformat()

if __name__ == '__main__':
    create_simple_test_case()
    print("Test data files created successfully!")
```

This comprehensive testing framework ensures the reliability, performance, and scientific validity of the Canopy-App model through multiple levels of automated testing.
