# Installation Guide

This guide will help you install and set up the Canopy-App modeling system.

## Prerequisites

### System Requirements

- **Operating System**: Linux, macOS
- **Compiler**: Modern Fortran compiler (gfortran, ifort, or similar)
- **Memory**: Minimum 4 GB RAM, 8 GB+ recommended
- **Storage**: At least 1 GB free disk space

### Required Dependencies

#### Fortran Compiler
```bash
# On Ubuntu/Debian
sudo apt-get install gfortran

# On CentOS/RHEL
sudo yum install gcc-gfortran

# On macOS with Homebrew
brew install gcc
```

#### NetCDF Library (Required for NetCDF I/O)
Canopy-App requires NetCDF-Fortran Libraries (`-lnetcdf -lnetcdff`) when using the 1D/2D NetCDF I/O Option (`infmt_opt=0`).

<!-- ```bash
# On Ubuntu/Debian
sudo apt-get install libnetcdf-dev libnetcdff-dev

# On CentOS/RHEL
sudo yum install netcdf-devel netcdf-fortran-devel

# On macOS with Homebrew
brew install netcdf netcdf-fortran

# On GMU Hopper (example module environment)
module load netcdf-c/4.7.4-vh netcdf-fortran/4.5.3-ff
``` -->

## Installation Methods

### Method 1: Clone from GitHub (Recommended)

```bash
# Clone the repository
git clone https://github.com/noaa-oar-arl/canopy-app.git
cd canopy-app

# Build with default settings (gfortran, NetCDF enabled)
make -C src
```

### Method 2: Download Release Archive

```bash
# Download and extract latest release
wget https://github.com/noaa-oar-arl/canopy-app/archive/main.zip
unzip main.zip
cd canopy-app-main

# Build
make -C src
```

## Build Configuration

Compilation options can be controlled with environment variables:

### Compiler Selection
- `FC=gfortran` (default) - GNU Fortran compiler
- `FC=ifort` - Intel Fortran compiler
- `FC=gfortran-11` - Specific GNU Fortran version
- `FC=/usr/bin/gfortran-11` - Full path to compiler

### Debug Options
- `DEBUG=0` (default) - No debug flags, optimized build
- `DEBUG=1` - Basic debug flags enabled
- `DEBUG=2` - Extensive debug flags including FPE traps and traceback

### NetCDF Options
- `NC=1` (default) - NetCDF support enabled
- `NC=0` - NetCDF support disabled (text I/O only)

## Build Examples

### Standard Build
```bash
# Default build (gfortran, optimized, NetCDF enabled)
make -C src
```

### Development Build with Debug
```bash
# Debug build with gfortran
DEBUG=1 NC=1 make -C src

# If FC is already set in environment, explicitly use gfortran
DEBUG=1 NC=1 FC=gfortran make -C src
```

### Intel Fortran Build
```bash
# Production build with Intel Fortran
FC=ifort make -C src

# Debug build with Intel Fortran
DEBUG=1 NC=1 FC=ifort make -C src
```

### Text-Only Build (No NetCDF)
```bash
# Build without NetCDF support
NC=0 make -C src
```

## Verification

After successful compilation, verify the installation:

```bash
# Check if executable was created
ls -la canopy

# # Check executable permissions
# ./canopy --help  # (if help option is implemented)

# # Or run with default settings
# ./canopy
# ```

## Troubleshooting

### Common Issues

#### NetCDF Not Found
```bash
# Check if nf-config is available
which nf-config
nf-config --fflags
nf-config --flibs

# If not found, install NetCDF development packages
# or set proper environment variables
```

#### Compiler Issues
```bash
# Check compiler version
gfortran --version
ifort --version

# Verify compiler can build simple Fortran programs
echo 'program test; print *, "Hello"; end program' > test.f90
gfortran test.f90 -o test
./test
rm test test.f90
```

#### Permission Issues
```bash
# Make sure you have write permissions
chmod +x canopy

# If building in a shared directory, check permissions
ls -la src/
```

### Environment Module Systems

For HPC systems using environment modules:

```bash
# Example for systems with NetCDF modules
module load netcdf-c/4.7.4
module load netcdf-fortran/4.5.3
module load compiler/gcc/9.3.0

# Build with loaded modules
make -C src
```

### Advanced Build Options

#### Custom Makefile Variables
```bash
# Override specific compiler flags
FCFLAGS="-O3 -march=native" make -C src

# Use custom NetCDF paths
NETCDF_INCDIR=/opt/netcdf/include \
NETCDF_LIBDIR=/opt/netcdf/lib \
make -C src
```

#### Parallel Build
```bash
# Use multiple cores for compilation
make -j4 -C src
```

## CMake Build (Recommended)

The new CMake-based build system supports cross-platform builds and advanced configuration. Use the provided `build.sh` script for a streamlined experience.

### Method 1: Build with `build.sh` (Recommended)

```bash
# From the project root
git clone https://github.com/noaa-oar-arl/canopy-app.git
cd canopy-app

# Build with default settings (gfortran, NetCDF enabled)
./build.sh
```

#### Common `build.sh` Options

- `--clean`         : Clean the build directory
- `--install`       : Run `make install` after build (installs to `install/` by default)
- `--no-modules`    : Skip environment module setup (for local/macOS/Linux builds)
- `-c <CMake opt>`  : Pass additional CMake options (can be used multiple times)
- `-t <target>`     : Specify build target (e.g., hera, macos, linux)

**Examples:**
```bash
# Clean and rebuild
./build.sh --clean && ./build.sh

# Build with NetCDF disabled
./build.sh -c "-DUSE_NETCDF=OFF"

# Build with debug flags
./build.sh -c "-DCANOPY_DEBUG_LEVEL=1"

# Install after build
./build.sh --install
```

### Method 2: Manual CMake Build

```bash
# From the project root
mkdir build
cd build
cmake ..
make -j4
make install  # Optional: install to ../install by default
```

## CMake Build Options

The CMake-based build system supports flexible configuration, including debug/release modes and NetCDF support.

### Common Build Customizations

- **Debug/Release Mode:**
  - `-DCANOPY_DEBUG_LEVEL=0` (Release, optimized, default)
  - `-DCANOPY_DEBUG_LEVEL=1` (Basic debug flags)
  - `-DCANOPY_DEBUG_LEVEL=2` (Extensive debug flags, FPE traps, traceback)
- **NetCDF Support:**
  - `-DUSE_NETCDF=ON` (default, enables NetCDF I/O)
  - `-DUSE_NETCDF=OFF` (disables NetCDF, text I/O only)

You can pass these options to `build.sh` using `-c` or directly to CMake:

**Examples:**

```bash
# Debug build with NetCDF enabled
./build.sh -c "-DCANOPY_DEBUG_LEVEL=1"

# Release build with NetCDF disabled
./build.sh -c "-DCANOPY_DEBUG_LEVEL=0" -c "-DUSE_NETCDF=OFF"

# Manual CMake build, debug, no NetCDF
mkdir build && cd build
cmake -DCANOPY_DEBUG_LEVEL=2 -DUSE_NETCDF=OFF ..
make -j4
```

You can combine these options as needed for your development or production workflow.

## Legacy Makefile Build (Deprecated)

> **Note:** The Makefile-based build is deprecated. Please use the new CMake-based build system and `./build.sh` for all new installations and development. The following instructions are for legacy users only and may be removed in a future release.

```bash
# Default build (gfortran, optimized, NetCDF enabled)
make -C src

# Debug build with gfortran
DEBUG=1 NC=1 make -C src

# If FC is already set in environment, explicitly use gfortran
DEBUG=1 NC=1 FC=gfortran make -C src

# Production build with Intel Fortran
FC=ifort make -C src

# Debug build with Intel Fortran
DEBUG=1 NC=1 FC=ifort make -C src

# Build without NetCDF support
NC=0 make -C src
```

## Verification

### Test Installation

Run the basic test to verify your installation:

```bash
# Run with example data
./canopy_app

# # Check version information
# ./canopy_app --version
# ```

### Expected Output

A successful installation should produce output similar to:

```
=========================================
    CANOPY-APP ATMOSPHERIC MODEL
=========================================
Version: 1.0
Build Date: 2024-XX-XX
Compiler: GNU Fortran
NetCDF Support: Enabled
=========================================
```

## Troubleshooting

### Common Issues

!!! warning "Compilation Errors"
    **Issue**: `gfortran: command not found`
    **Solution**: Install a Fortran compiler (see prerequisites)

!!! warning "NetCDF Errors"
    **Issue**: `fatal error: netcdf.mod: No such file or directory`
    **Solution**: Install NetCDF development libraries or compile without NetCDF

!!! warning "Permission Errors"
    **Issue**: `Permission denied` when running executable
    **Solution**: Make the file executable: `chmod +x canopy_app`

### Environment Variables

You may need to set environment variables for libraries:

```bash
# For NetCDF
export NETCDF_ROOT=/usr/local
export LD_LIBRARY_PATH=$NETCDF_ROOT/lib:$LD_LIBRARY_PATH

# For Intel Fortran
export INTEL_LICENSE_FILE=/opt/intel/licenses
```

### Getting Help

If you encounter issues during installation:

1. Check the [troubleshooting section](../user-guide/troubleshooting.md)
2. Search [existing issues](https://github.com/canopy-app/canopy-app/issues)
3. Create a [new issue](https://github.com/canopy-app/canopy-app/issues/new) with:
   - Your operating system
   - Compiler version
   - Complete error messages
   - Installation method used

## Next Steps

Once installation is complete:

1. 📖 Read the [Quick Start Guide](quickstart.md)
2. ⚙️ Configure your [model settings](configuration.md)
3. 🚀 Run your [first simulation](../examples/basic.md)
