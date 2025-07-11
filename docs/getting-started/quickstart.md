# Quick Start Guide

Get up and running with Canopy-App in just a few steps!

## Prerequisites Check

Before starting, ensure you have:
- [x] Fortran compiler (gfortran recommended)
- [x] NetCDF library (optional but recommended)
- [x] Input data files

## 1. Clone and Compile

```bash
# Clone the repository
git clone https://github.com/canopy-app/canopy-app.git
cd canopy-app

# Compile the model
cd src
make

# Check compilation was successful
ls -la canopy
```

## 2. Prepare Input Files

The model requires several input files:

### Required Files
- `namelist.canopy` - Main configuration file
- Meteorological input data (NetCDF or text format)
- Point file data (optional)

### Example Setup
```bash
# Copy example input files
cp input/namelist.canopy .
cp input/gfs.t12z.20220701.sfcf000.canopy.nc .
```

## 3. Configure the Model

Edit the namelist file to match your setup:

```fortran
&canopy_inputs
    file_vars = 'gfs.t12z.20220701.sfcf000.canopy.nc'
    file_out = 'output/canopy_output'
/
```

## 4. Run the Model

```bash
# Run with modified namelist file
./canopy_app



## 5. Check Output

```bash
# List output files
ls -la *.nc *.txt

# Quick check of NetCDF output
ncdump -h canopy_output.nc
```

## Example Output

A successful run will produce:

```
 CANOPY-APP: Starting model execution
 Reading namelist: namelist.canopy
 Processing date: 20220701, time: 00
 Computing canopy meteorology...
 Computing radiation transfer...
 Computing biogenic emissions...
 Computing dry deposition...
 Writing output files...
 CANOPY-APP: Model execution completed successfully
```

## Next Steps

- 📖 [Learn about configuration options](configuration.md)
- 🔍 [Explore the User Guide](../user-guide/overview.md)
- 📊 [Try more examples](../examples/basic.md)
- 🐛 [Troubleshooting common issues](../user-guide/troubleshooting.md)

## Common Issues

!!! warning "Compilation Errors"
    If you encounter compilation errors, check:
    - Fortran compiler version (`gfortran --version`)
    - NetCDF library installation
    - Makefile paths and flags

!!! tip "Performance Tips"
    - Use optimized compiler flags for production runs
    - Consider parallel processing for large domains
    - Monitor memory usage for extensive simulations
