# Input Files

The Canopy-App model requires several input files to run successfully. This section describes the format and content of each required input file.

## Namelist Configuration File

The primary configuration file is `namelist.canopy`, which contains all model parameters and settings.

### File Format

The namelist file uses Fortran namelist format with the following structure:

```fortran
&CANOPY_OPTIONS
 file_vars = 'input/gfs.t12z.20220630.sfcf023.canopy.nc'
 infmt_opt = 1
 nlat = 1
 nlon = 1
 ntime = 1
 time_start = '2022-07-01-11:00:00.0000'

/
```

### Key Parameters

| Parameter | Description | Units | Default |
|-----------|-------------|-------|---------|
| `infmt_opt` | Input format option (0=netCDF, 1=text) | - | 1 |
| `nlat` | Number of latitude points | - | 1 |
| `nlon` | Number of longitude points | - | 1 |
| `ntime` | Number of time steps | - | 1 |
| `time_start` | Start time (YYYY-MM-DD_HH:MM:SS) | - | - |


## Meteorological Input Files

### NetCDF Format

When using `infmt_opt = 0`, the model expects netCDF files containing meteorological variables:

- **File naming**: `gfs.tXXz.YYYYMMDD.sfcfXXX.canopy.nc`
- **Example required variables**:
  - `tmp2m`: 2-meter temperature (K)
  - `spfh2m`: 2-meter specific humidity (kg/kg)
  - `pressfc`: Surface pressure (Pa)
  - `fricv`: Friction velocity (m/s)
  - `ugrd10m`: U-component of 10-meter wind speed (m/s)
  - `vgrd10m`: V-component of 10-meter wind speed (m/s)
  - `dswrf`: Instantaneous downward shortwave radiation at surface (W/m²)

### Text Format

When using `infmt_opt = 1`, the model reads text files with the following format:

```
# Time: 2022-06-30_12:00:00
# Lat: 40.0, Lon: -80.0
lat      lon      ch          ugrd10m  vgrd10m
34.0   272.1   20.9      -0.75          0.33
```

## Vegetation Parameters

Vegetation Parameters are assumed to come through input or configuration at the moment

<!-- ## File Validation

The model performs input validation checks:

1. **Format validation**: Ensures files match expected format
2. **Range validation**: Checks that values are within physical limits
3. **Consistency validation**: Verifies compatibility between different inputs
4. **Time validation**: Ensures time stamps are consistent -->

## Troubleshooting Input Files

Common issues and solutions:

- **Missing variables**: Ensure all required variables are present in netCDF files
- **Incorrect dimensions**: Check that spatial and temporal dimensions match namelist settings
- **Invalid values**: Look for NaN, infinity, or out-of-range values
- **Time format errors**: Use exact format YYYY-MM-DD_HH:MM:SS

For more troubleshooting help, see the [Troubleshooting Guide](troubleshooting.md).
