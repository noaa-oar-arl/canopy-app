# Performance Optimization

This guide provides comprehensive information on optimizing the performance of the Canopy-App model for various computational environments and use cases.

## Performance Analysis

### Profiling Tools

#### Fortran Profiling

Use `gprof` for detailed performance analysis:

```bash
# Compile with profiling enabled
make clean
make FFLAGS="-pg -O2"

# Run the model
./canopy_app.exe

# Generate profile report
gprof canopy_app.exe gmon.out > profile.txt

# Analyze the results
head -30 profile.txt
```

#### Intel VTune (for Intel Compiler)

```bash
# Compile with debug symbols
make FFLAGS="-g -O3"

# Profile with VTune
vtune -collect hotspots -result-dir vtune_results ./canopy_app.exe

# View results
vtune-gui vtune_results
```

#### System-Level Monitoring

```bash
# Monitor CPU and memory usage
top -p $(pgrep canopy_app)

# Monitor I/O performance
iostat -x 1

# Monitor memory bandwidth
perf stat -e cache-misses,cache-references ./canopy_app.exe
```

### Performance Bottlenecks

Common performance bottlenecks identified through profiling:

1. **I/O Operations** (typically 20-40% of runtime)
2. **Radiation Calculations** (15-25% of runtime)
3. **Photosynthesis Calculations** (10-20% of runtime)
4. **Memory Allocation/Deallocation** (5-15% of runtime)
5. **Turbulence Calculations** (5-10% of runtime)

## Compiler Optimization

### GCC/Gfortran Optimization

```makefile
# Basic optimization
FFLAGS = -O3 -march=native -mtune=native

# Advanced optimization
FFLAGS = -O3 -march=native -mtune=native \
         -funroll-loops -ffast-math \
         -flto -fwhole-program

# Debugging-friendly optimization
FFLAGS = -O2 -g -march=native -mtune=native

# Profile-guided optimization
FFLAGS_PGO_GEN = -O3 -fprofile-generate
FFLAGS_PGO_USE = -O3 -fprofile-use
```

### Intel Fortran Optimization

```makefile
# Basic optimization
FFLAGS = -O3 -xHost -ipo

# Advanced optimization
FFLAGS = -O3 -xHost -ipo -fast \
         -no-prec-div -fp-model fast=2

# Vectorization reports
FFLAGS = -O3 -vec-report=2 -opt-report=2

# Profile-guided optimization
FFLAGS_PGO_GEN = -O3 -prof-gen
FFLAGS_PGO_USE = -O3 -prof-use
```

### Architecture-Specific Optimization

```bash
# Check CPU capabilities
cat /proc/cpuinfo | grep flags

# Optimize for specific CPU
# Intel Skylake
FFLAGS="-O3 -march=skylake -mtune=skylake"

# AMD Zen2
FFLAGS="-O3 -march=znver2 -mtune=znver2"

# Generic modern CPU
FFLAGS="-O3 -march=native -mtune=native"
```

## Memory Optimization

### Memory Layout Optimization

#### Array Access Patterns

Optimize for cache efficiency:

```fortran
! BAD: Non-contiguous memory access
do k = 1, nlev
  do i = 1, nlon
    do j = 1, nlat
      temp(i,j,k) = temp(i,j,k) + heating_rate
    end do
  end do
end do

! GOOD: Contiguous memory access
do k = 1, nlev
  do j = 1, nlat
    do i = 1, nlon
      temp(i,j,k) = temp(i,j,k) + heating_rate
    end do
  end do
end do
```

#### Memory Alignment

```fortran
! Align arrays to cache boundaries
real(r8), allocatable :: temp(:,:,:)
!dir$ attributes align:64 :: temp

! Use contiguous array attribute
real(r8), contiguous, pointer :: temp_ptr(:,:,:)
```

### Memory Pool Management

Implement memory pools for frequently allocated/deallocated arrays:

```fortran
module memory_pool_mod
  implicit none
  private

  type :: memory_pool_type
    real(r8), allocatable :: pool(:)
    integer :: pool_size
    integer :: current_pos
  end type

  public :: memory_pool_type, init_pool, get_memory, return_memory

contains

  subroutine init_pool(pool, size_mb)
    type(memory_pool_type), intent(out) :: pool
    integer, intent(in) :: size_mb

    pool%pool_size = size_mb * 1024 * 1024 / 8  ! Convert MB to real(r8) elements
    allocate(pool%pool(pool%pool_size))
    pool%current_pos = 1
  end subroutine

end module
```

### Large Memory Optimizations

For systems with limited memory:

```fortran
! Process data in chunks
integer, parameter :: chunk_size = 1000

do chunk_start = 1, total_points, chunk_size
  chunk_end = min(chunk_start + chunk_size - 1, total_points)

  ! Process chunk
  call process_chunk(data(chunk_start:chunk_end))
end do
```

## Parallel Performance

### OpenMP Optimization

#### Thread Scaling

```bash
# Test different thread counts
for threads in 1 2 4 8 16; do
  export OMP_NUM_THREADS=$threads
  echo "Testing with $threads threads:"
  time ./canopy_app.exe
done
```

#### Load Balancing

```fortran
! Use dynamic scheduling for uneven work
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,10) PRIVATE(i,j,k)
do ipoint = 1, total_points
  call expensive_calculation(ipoint)
end do
!$OMP END PARALLEL DO

! Use guided scheduling for decreasing work
!$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,j,k)
do k = 1, nlev
  do j = 1, nlat
    do i = 1, nlon
      call point_calculation(i,j,k)
    end do
  end do
end do
!$OMP END PARALLEL DO
```

#### Thread Affinity

```bash
# Bind threads to cores
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Explicit thread affinity
export GOMP_CPU_AFFINITY="0-7"
```

#### Reducing False Sharing

```fortran
! BAD: False sharing
real(r8) :: thread_sum(max_threads)
!$OMP PARALLEL
thread_id = omp_get_thread_num()
thread_sum(thread_id) = calculate_sum()
!$OMP END PARALLEL

! GOOD: Avoid false sharing with padding
type :: padded_sum_type
  real(r8) :: sum
  real(r8) :: padding(15)  ! Pad to cache line size
end type
type(padded_sum_type) :: thread_sums(max_threads)
```

### NUMA Optimization

For Non-Uniform Memory Access systems:

```bash
# Check NUMA topology
numactl --hardware

# Run with NUMA awareness
numactl --cpunodebind=0 --membind=0 ./canopy_app.exe

# Interleave memory across nodes
numactl --interleave=all ./canopy_app.exe
```

## I/O Optimization

### NetCDF Optimization

#### Chunking Strategy

```fortran
! Optimize chunking for access patterns
integer :: chunk_sizes(4) = [32, 32, 10, 1]  ! [lon, lat, lev, time]

status = nf90_def_var_chunking(ncid, varid, NF90_CHUNKED, chunk_sizes)
```

#### Compression Settings

```fortran
! Enable compression
status = nf90_def_var_deflate(ncid, varid, shuffle=1, deflate=1, deflate_level=6)

! Use appropriate fill values
status = nf90_def_var_fill(ncid, varid, no_fill=0, fill_value=missing_value)
```

#### Parallel I/O

```fortran
! Enable parallel I/O (if available)
status = nf90_create_par(filename, ior(NF90_CLOBBER, NF90_MPIIO), &
                        MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
```

### File System Optimization

#### Temporary Files

```bash
# Use fast local storage for temporary files
export TMPDIR=/tmp/fast_storage

# Or use memory filesystem
export TMPDIR=/dev/shm
```

#### I/O Patterns

```fortran
! Write large blocks rather than small frequent writes
integer, parameter :: buffer_size = 1000000
real(r8) :: buffer(buffer_size)
integer :: buffer_pos = 0

do i = 1, n_data_points
  buffer_pos = buffer_pos + 1
  buffer(buffer_pos) = data(i)

  if (buffer_pos == buffer_size) then
    call write_buffer(buffer, buffer_pos)
    buffer_pos = 0
  end if
end do

! Write remaining data
if (buffer_pos > 0) then
  call write_buffer(buffer, buffer_pos)
end if
```

## Algorithm Optimization

### Numerical Methods

#### Fast Mathematical Functions

```fortran
! Use lookup tables for expensive functions
module fast_math_mod
  implicit none
  private

  integer, parameter :: lut_size = 10000
  real(r8) :: exp_lut(lut_size)
  real(r8) :: log_lut(lut_size)

  public :: fast_exp, fast_log

contains

  function fast_exp(x) result(result)
    real(r8), intent(in) :: x
    real(r8) :: result
    integer :: index

    ! Clamp and scale input
    index = max(1, min(lut_size, int((x + 10.0) * lut_size / 20.0)))
    result = exp_lut(index)
  end function

end module
```

#### Vectorization-Friendly Code

```fortran
! GOOD: Vectorizable loop
do i = 1, n
  result(i) = sqrt(a(i)*a(i) + b(i)*b(i))
end do

! BAD: Non-vectorizable due to dependency
do i = 2, n
  result(i) = result(i-1) + data(i)
end do
```

### Physics Optimizations

#### Radiation Calculations

```fortran
! Pre-compute expensive trigonometric functions
module radiation_opt_mod
  real(r8), save :: cos_zenith
  real(r8), save :: sin_zenith
  logical, save :: trig_computed = .false.

contains

  subroutine compute_solar_angles(zenith_angle)
    real(r8), intent(in) :: zenith_angle

    if (.not. trig_computed) then
      cos_zenith = cos(zenith_angle)
      sin_zenith = sin(zenith_angle)
      trig_computed = .true.
    end if
  end subroutine

end module
```

#### Photosynthesis Calculations

```fortran
! Cache temperature-dependent parameters
type :: photosynthesis_cache_type
  real(r8) :: last_temp = -999.0
  real(r8) :: vcmax_temp
  real(r8) :: jmax_temp
  real(r8) :: kc_temp
  real(r8) :: ko_temp
end type

subroutine photosynthesis_optimized(temp, ppfd, co2, result, cache)
  real(r8), intent(in) :: temp, ppfd, co2
  real(r8), intent(out) :: result
  type(photosynthesis_cache_type), intent(inout) :: cache

  ! Only recalculate if temperature changed significantly
  if (abs(temp - cache%last_temp) > 0.1) then
    call update_temperature_params(temp, cache)
    cache%last_temp = temp
  end if

  ! Use cached parameters
  result = calculate_photosynthesis(ppfd, co2, cache)
end subroutine
```

## Performance Tuning Guidelines

### Development Phase

1. **Profile early and often** - Use profiling tools throughout development
2. **Optimize hot spots first** - Focus on code that consumes most time
3. **Measure before and after** - Quantify performance improvements
4. **Consider algorithmic improvements** - Often more effective than micro-optimizations

### Compilation Phase

1. **Use appropriate optimization levels** - Balance performance vs. compile time
2. **Enable processor-specific optimizations** - Use `-march=native` or specific targets
3. **Consider profile-guided optimization** - For production builds
4. **Test different compilers** - Intel, GCC, and PGI may perform differently

### Runtime Phase

1. **Choose optimal thread count** - Usually equals physical CPU cores
2. **Monitor memory usage** - Avoid swapping at all costs
3. **Use appropriate data formats** - Binary formats are faster than text
4. **Consider batch processing** - For multiple simulations

## Performance Monitoring

### Automated Benchmarking

```bash
#!/bin/bash
# performance_test.sh

# Test configurations
CONFIGS=("config1.nml" "config2.nml" "config3.nml")
THREADS=(1 2 4 8 16)

for config in "${CONFIGS[@]}"; do
  for threads in "${THREADS[@]}"; do
    echo "Testing $config with $threads threads"
    export OMP_NUM_THREADS=$threads

    # Run 3 times and take average
    times=()
    for run in {1..3}; do
      start_time=$(date +%s.%N)
      ./canopy_app.exe $config > /dev/null 2>&1
      end_time=$(date +%s.%N)
      runtime=$(echo "$end_time - $start_time" | bc)
      times+=($runtime)
    done

    # Calculate average
    avg_time=$(echo "${times[@]}" | awk '{sum=0; for(i=1;i<=NF;i++)sum+=$i; print sum/NF}')
    echo "Average runtime: $avg_time seconds"
  done
done
```

### Performance Regression Testing

```python
#!/usr/bin/env python3
# performance_regression.py

import subprocess
import time
import json
import sys

def run_benchmark(config_file, num_threads=4):
    """Run benchmark and return execution time"""
    env = {'OMP_NUM_THREADS': str(num_threads)}

    start_time = time.time()
    result = subprocess.run(['./canopy_app.exe', config_file],
                          env=env, capture_output=True)
    end_time = time.time()

    if result.returncode != 0:
        raise RuntimeError(f"Benchmark failed: {result.stderr}")

    return end_time - start_time

def main():
    baseline_file = 'performance_baseline.json'

    # Run current benchmark
    current_time = run_benchmark('benchmark_config.nml')

    # Load baseline
    try:
        with open(baseline_file, 'r') as f:
            baseline = json.load(f)
        baseline_time = baseline['execution_time']
    except FileNotFoundError:
        # Create baseline
        baseline = {'execution_time': current_time}
        with open(baseline_file, 'w') as f:
            json.dump(baseline, f)
        print(f"Created baseline: {current_time:.2f} seconds")
        return

    # Check for regression
    regression_threshold = 1.1  # 10% slower
    if current_time > baseline_time * regression_threshold:
        print(f"PERFORMANCE REGRESSION DETECTED!")
        print(f"Baseline: {baseline_time:.2f}s, Current: {current_time:.2f}s")
        print(f"Slowdown: {(current_time/baseline_time-1)*100:.1f}%")
        sys.exit(1)

    print(f"Performance OK: {current_time:.2f}s (baseline: {baseline_time:.2f}s)")

if __name__ == '__main__':
    main()
```

## Hardware-Specific Optimizations

### Intel Xeon Processors

```makefile
# Intel Xeon optimization
FFLAGS = -O3 -xCORE-AVX512 -qopt-zmm-usage=high -ipo
```

### AMD EPYC Processors

```makefile
# AMD EPYC optimization
FFLAGS = -O3 -march=znver2 -mtune=znver2 -mavx2
```

### ARM Processors

```makefile
# ARM Neoverse optimization
FFLAGS = -O3 -mcpu=neoverse-n1 -mtune=neoverse-n1
```

By following these optimization guidelines, you can significantly improve the performance of the Canopy-App model for your specific use case and hardware configuration.
