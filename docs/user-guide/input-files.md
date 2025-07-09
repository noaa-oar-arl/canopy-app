# Input Files

The Canopy-App model requires several input files to run successfully. This section describes the format and content of each required input file.

## Namelist Configuration File

The primary configuration file is `namelist.canopy`, which contains all model parameters and settings.

### File Format

The namelist file uses Fortran namelist format with the following structure:

```fortran
&CANOPY_OPTIONS
 file_vars = 'namelist.canopy'
 infmt_opt = 1
 nlat = 1
 nlon = 1
 ntime = 1
 time_start = '2022-06-30_12:00:00'
 dx = 100.0
 dy = 100.0
 dz_top = 3.0
/
```

### Key Parameters

| Parameter | Description | Units | Default |
|-----------|-------------|-------|---------|
| `infmt_opt` | Input format option (1=netCDF, 2=text) | - | 1 |
| `nlat` | Number of latitude points | - | 1 |
| `nlon` | Number of longitude points | - | 1 |
| `ntime` | Number of time steps | - | 1 |
| `time_start` | Start time (YYYY-MM-DD_HH:MM:SS) | - | - |
| `dx` | Grid spacing in x-direction | m | 100.0 |
| `dy` | Grid spacing in y-direction | m | 100.0 |
| `dz_top` | Top layer thickness | m | 3.0 |

## Meteorological Input Files

### NetCDF Format

When using `infmt_opt = 1`, the model expects netCDF files containing meteorological variables:

- **File naming**: `gfs.tXXz.YYYYMMDD.sfcfXXX.canopy.nc`
- **Required variables**:
  - `TEMP_2M`: 2-meter temperature (K)
  - `QV_2M`: 2-meter specific humidity (kg/kg)
  - `PRES_SFC`: Surface pressure (Pa)
  - `USTAR`: Friction velocity (m/s)
  - `WSPD_10M`: 10-meter wind speed (m/s)
  - `WDIR_10M`: 10-meter wind direction (degrees)
  - `SRAD_TOA`: Top-of-atmosphere solar radiation (W/m²)

### Text Format

When using `infmt_opt = 2`, the model reads text files with the following format:

```
# Time: 2022-06-30_12:00:00
# Lat: 40.0, Lon: -80.0
TEMP_2M    QV_2M      PRES_SFC   USTAR     WSPD_10M  WDIR_10M  SRAD_TOA
295.15     0.012      101325.0   0.45      5.2       270.0     850.0
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
