# Troubleshooting

This guide helps you diagnose and resolve common issues when running the Canopy-App model.

## Common Error Messages

### Compilation Errors

#### Missing Dependencies

```
Error: netcdf.mod not found
```

**Solution:**
```bash
# Install NetCDF Fortran development libraries
# Ubuntu/Debian:
sudo apt-get install libnetcdff-dev

# CentOS/RHEL:
sudo yum install netcdf-fortran-devel

# macOS with Homebrew:
brew install netcdf
```

#### Compiler Version Issues

```
Error: Fortran 2008 features not supported
```

**Solution:**
- Ensure you're using gfortran 4.8+ or ifort 15.0+
- Update compiler flags in Makefile:
  ```makefile
  FFLAGS = -std=f2008 -fall-intrinsics
  ```

### Runtime Errors

#### Memory Allocation Errors

```
Error: Out of memory during array allocation
```

**Causes and Solutions:**

1. **Grid too large for available memory**
   ```bash
   # Check available memory
   free -h
   # Reduce grid size in namelist
   nlat = 50  # instead of 500
   nlon = 50  # instead of 500
   ```

2. **Memory leaks**
   ```bash
   # Compile with debugging
   make clean
   make FFLAGS="-g -fbounds-check -fcheck=all"
   ```

3. **System memory limits**
   ```bash
   # Check and increase limits
   ulimit -a
   ulimit -s unlimited  # stack size
   ulimit -v unlimited  # virtual memory
   ```

#### File I/O Errors

```
Error: Cannot open file 'input.nc'
```

**Solutions:**

1. **Check file existence and permissions**
   ```bash
   ls -la input.nc
   chmod 644 input.nc
   ```

2. **Verify file format**
   ```bash
   ncdump -h input.nc  # Check netCDF structure
   file input.nc       # Check file type
   ```

3. **Check file path in namelist**
   ```fortran
   &CANOPY_OPTIONS
    file_vars = './input/namelist.canopy'  ! Use relative path
   /
   ```

#### Numerical Instability

```
Warning: Solution did not converge
Error: NaN values detected in output
```

**Solutions:**

1. **Check input data quality**
   ```python
   import netCDF4 as nc
   ds = nc.Dataset('input.nc')
   temp = ds.variables['temp'][:]
   print(f"Temperature range: {temp.min():.2f} to {temp.max():.2f} K")
   print(f"NaN count: {np.isnan(temp).sum()}")
   ```


## Performance Issues

### Slow Execution

#### Profiling

```bash
# Compile with profiling
make FFLAGS="-pg -O2"

# Run simulation
./canopy_app.exe

# Analyze profile
gprof canopy_app.exe gmon.out > profile.txt
```

#### Common Bottlenecks

1. **I/O Performance**
   ```fortran
   ! Use netCDF instead of text files
   infmt_opt = 0  ! netCDF input
   ```

2. **Memory Access Patterns**
   ```fortran
   ! Optimize loop order for cache efficiency
   do k = 1, nlev
     do j = 1, nlat
       do i = 1, nlon
         ! Process data(i,j,k)
       end do
     end do
   end do
   ```

3. **Compiler Optimization**
   ```bash
   make FFLAGS="-O3 -march=native -fopenmp"
   ```

### Memory Usage

#### Memory Monitoring

```bash
# Monitor memory during execution
top -p $(pgrep canopy_app)

# Or use memory profiling
valgrind --tool=massif ./canopy_app.exe
```

#### Memory Optimization

1. **Reduce precision where appropriate**
   ```fortran
   real(kind=real32) :: temp_array  ! Single precision
   ```

2. **Use dynamic allocation**
   ```fortran
   allocate(temp_array(nlat, nlon))
   ! ... use array ...
   deallocate(temp_array)
   ```

## Input Data Issues

### Data Validation

#### Check Variable Ranges

```python
import netCDF4 as nc
import numpy as np

def validate_input(filename):
    ds = nc.Dataset(filename, 'r')

    # Check temperature
    temp = ds.variables['temp'][:]
    if temp.min() < 150 or temp.max() > 350:
        print(f"Warning: Temperature out of range: {temp.min():.1f} - {temp.max():.1f} K")

    # Check humidity
    qv = ds.variables['qv'][:]
    if qv.min() < 0 or qv.max() > 0.1:
        print(f"Warning: Humidity out of range: {qv.min():.6f} - {qv.max():.6f} kg/kg")

    # Check for missing values
    for var in ['temp', 'qv', 'u', 'v']:
        data = ds.variables[var][:]
        nan_count = np.isnan(data).sum()
        if nan_count > 0:
            print(f"Warning: {nan_count} NaN values in {var}")

    ds.close()

validate_input('input.nc')
```

#### Time Coordinate Issues

```python
def check_time_coords(filename):
    ds = nc.Dataset(filename, 'r')
    time = ds.variables['time']

    print(f"Time units: {time.units}")
    print(f"Time range: {time[:].min()} to {time[:].max()}")

    # Check for regular spacing
    dt = np.diff(time[:])
    if not np.allclose(dt, dt[0]):
        print("Warning: Irregular time spacing detected")

    ds.close()
```

### Data Preprocessing

#### Fill Missing Values

```python
import numpy as np
from scipy.interpolate import griddata

def fill_missing_values(data):
    """Fill NaN values using spatial interpolation"""
    if np.isnan(data).any():
        mask = ~np.isnan(data)
        coords = np.array(np.where(mask)).T
        values = data[mask]
        missing_coords = np.array(np.where(~mask)).T

        if len(missing_coords) > 0:
            filled = griddata(coords, values, missing_coords, method='linear')
            data[~mask] = filled

    return data
```

#### Smooth Noisy Data

```python
from scipy.ndimage import gaussian_filter

def smooth_data(data, sigma=1.0):
    """Apply Gaussian smoothing to reduce noise"""
    return gaussian_filter(data, sigma=sigma)
```

## Model Configuration Issues

### Namelist Validation

#### Check Parameter Ranges

```python
def validate_namelist(namelist_file):
    """Validate namelist parameters"""
    with open(namelist_file, 'r') as f:
        content = f.read()

    # Extract parameters (simplified)
    import re

    # Check grid dimensions
    nlat_match = re.search(r'nlat\s*=\s*(\d+)', content)
    nlon_match = re.search(r'nlon\s*=\s*(\d+)', content)

    if nlat_match and nlon_match:
        nlat, nlon = int(nlat_match.group(1)), int(nlon_match.group(1))
        if nlat * nlon > 10000:
            print(f"Warning: Large grid size ({nlat}x{nlon}) may cause memory issues")

    # Check time step
    dt_match = re.search(r'dt\s*=\s*([\d.]+)', content)
    if dt_match:
        dt = float(dt_match.group(1))
        if dt > 3600:
            print(f"Warning: Large time step ({dt}s) may cause instability")

validate_namelist('namelist.canopy')
```

### Physics Parameter Issues

#### Unrealistic Vegetation Parameters

```fortran
! Check vegetation parameters in namelist
&VEGETATION
 lai = 4.0      ! Reasonable for forest (0.1-10.0)
 canht = 20.0   ! Reasonable height in meters (0.1-50.0)
 z0 = 2.0       ! Should be ~10% of canht
/
```

## Debugging Techniques

<!-- TODO: Fix this section
### Enable Debug Output

```fortran
&CANOPY_OPTIONS
 debug_level = 2  ! 0=none, 1=basic, 2=detailed, 3=verbose
/
``` -->

### Add Debug Prints

```fortran
subroutine my_subroutine(temp, qv)
  real, intent(in) :: temp, qv

  ! Add debug output
  if (debug_level >= 2) then
    write(*,*) 'DEBUG: my_subroutine called with temp=', temp, 'qv=', qv
  end if

  ! Check for invalid values
  if (temp < 0.0 .or. temp > 400.0) then
    write(*,*) 'ERROR: Invalid temperature:', temp
    stop
  end if
end subroutine
```

### Use Compiler Debug Flags

```bash
# Debug compilation
make clean
make FFLAGS="-g -O0 -fbounds-check -fcheck=all -Wall"

# Run with debugger
gdb ./canopy_app.exe
(gdb) run
# ... if crash occurs ...
(gdb) bt  # backtrace
(gdb) print variable_name
```

## Getting Help

### Log Files

Always check log files for error messages:

```bash
# Check standard output
cat canopy_output.log

# Check error output
cat canopy_error.log

# Check system logs
dmesg | grep -i "canopy\|memory\|segfault"
```

### System Information

Gather system information for bug reports:

```bash
# System info
uname -a
cat /proc/version

# Compiler version
gfortran --version
ifort --version

# Library versions
nc-config --version
nf-config --version

# Memory and disk space
free -h
df -h
```

### Reporting Issues

When reporting issues, include:

1. **Error message** (complete text)
2. **System information** (OS, compiler, libraries)
3. **Input files** (namelist and data files)
4. **Steps to reproduce**
5. **Expected vs. actual behavior**

### Community Resources

- **GitHub Issues**: Report bugs and feature requests 
- **Documentation**: Check latest online documentation

## Quick Reference

### Common Commands

```bash
# Check file format
file input.nc
ncdump -h input.nc

# Monitor resources
top -p $(pgrep canopy_app)
watch -n 1 'free -h'

# Debug compilation
make clean && make FFLAGS="-g -fbounds-check"

# Memory check
valgrind --leak-check=full ./canopy_app.exe
```

### Emergency Stops

```bash
# Kill running simulation
pkill -f canopy_app

# Force kill if necessary
pkill -9 -f canopy_app

# Check for zombie processes
ps aux | grep canopy_app
```
