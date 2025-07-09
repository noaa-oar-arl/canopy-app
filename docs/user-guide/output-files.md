# Output Files

The Canopy-App model generates various output files containing simulation results. This section describes the format and content of each output file type.

## Output File Types

### NetCDF Output Files

The model generates netCDF files with comprehensive simulation results:

<!-- #### Primary Output File: `canopy_output.nc`

Contains the main simulation variables:

**Dimensions:**
- `time`: Number of time steps
- `lev`: Number of vertical levels
- `lat`: Number of latitude points
- `lon`: Number of longitude points

**Variables:**

| Variable | Description | Units | Dimensions |
|----------|-------------|-------|------------|
| `time` | Time coordinate | hours since start | (time) |
| `lev` | Vertical level coordinate | m | (lev) |
| `lat` | Latitude coordinate | degrees_north | (lat) |
| `lon` | Longitude coordinate | degrees_east | (lon) |
| `temp` | Air temperature profile | K | (time,lev,lat,lon) |
| `qv` | Water vapor mixing ratio | kg/kg | (time,lev,lat,lon) |
| `u` | U-component wind speed | m/s | (time,lev,lat,lon) |
| `v` | V-component wind speed | m/s | (time,lev,lat,lon) |
| `tke` | Turbulent kinetic energy | m²/s² | (time,lev,lat,lon) |

#### Canopy-Specific Variables

| Variable | Description | Units | Dimensions |
|----------|-------------|-------|------------|
| `lai` | Leaf area index | m²/m² | (time,lat,lon) |
| `canht` | Canopy height | m | (time,lat,lon) |
| `ppfd_sun` | PPFD for sunlit leaves | μmol/m²/s | (time,lev,lat,lon) |
| `ppfd_shade` | PPFD for shaded leaves | μmol/m²/s | (time,lev,lat,lon) |
| `tleaf_sun` | Sunlit leaf temperature | K | (time,lev,lat,lon) |
| `tleaf_shade` | Shaded leaf temperature | K | (time,lev,lat,lon) |
| `gsw_sun` | Sunlit stomatal conductance | mol/m²/s | (time,lev,lat,lon) |
| `gsw_shade` | Shaded stomatal conductance | mol/m²/s | (time,lev,lat,lon) |

### Text Output Files

For simple analysis, the model can output text files:

#### Point Output: `point_file_YYYYMMDD.sfcfXXX.txt`

Contains time series data for single points:

```
# Canopy-App Point Output
# Time: 2022-06-30_12:00:00, Lat: 40.0, Lon: -80.0
# Level(m)  Temp(K)  QV(kg/kg)  U(m/s)  V(m/s)  TKE(m2/s2)
    0.5     295.15    0.0120     2.1     0.8      0.45
    1.5     294.98    0.0118     2.3     1.0      0.52
    2.5     294.82    0.0116     2.6     1.2      0.58
    ...
```

#### Canopy Variables: `point_file_canvars_YYYYMMDD.sfcfXXX.txt`

Contains canopy-specific variables:

```
# Canopy Variables Output
# Time: 2022-06-30_12:00:00, Lat: 40.0, Lon: -80.0
LAI:          4.50
CANHT:       15.20
PPFD_SUN:   1250.30
PPFD_SHADE:  185.70
TLEAF_SUN:   298.45
TLEAF_SHADE: 296.20
GSW_SUN:     0.185
GSW_SHADE:   0.045
``` -->

## Chemical Output (if enabled)

### Species Concentrations

When chemistry is enabled, additional files contain species concentrations:

#### `canopy_chem.nc`

Contains chemical species profiles:

| Variable | Description | Units | Dimensions |
|----------|-------------|-------|------------|
| `co2` | CO₂ concentration | ppm | (time,lev,lat,lon) |
| `h2o` | H₂O concentration | ppm | (time,lev,lat,lon) |
| `o3` | O₃ concentration | ppb | (time,lev,lat,lon) |
| `no` | NO concentration | ppb | (time,lev,lat,lon) |
| `no2` | NO₂ concentration | ppb | (time,lev,lat,lon) |
| `so2` | SO₂ concentration | ppb | (time,lev,lat,lon) |

### Emission Rates

#### `canopy_emis.nc`

Contains biogenic emission rates:

| Variable | Description | Units | Dimensions |
|----------|-------------|-------|------------|
| `isop_emis` | Isoprene emission rate | μg/m²/s | (time,lat,lon) |
| `mono_emis` | Monoterpene emission rate | μg/m²/s | (time,lat,lon) |
| `sesq_emis` | Sesquiterpene emission rate | μg/m²/s | (time,lat,lon) |
| `ovoc_emis` | Other VOC emission rate | μg/m²/s | (time,lat,lon) |

## Diagnostic Files

### Model Performance

#### `canopy_timing.txt`

Contains timing information:

```
Canopy-App Performance Report
=============================
Total Runtime:        00:02:35.123
Initialization:       00:00:05.234
Main Loop:           00:02:25.456
  - Meteorology:     00:00:45.123
  - Canopy Physics:  00:01:15.234
  - Chemistry:       00:00:20.456
  - I/O:            00:00:04.643
Finalization:        00:00:04.433

Memory Usage:
Peak Memory:         1.23 GB
Average Memory:      0.89 GB
```

#### `canopy_log.txt`

Contains detailed runtime information:

```
2022-06-30 12:00:00 [INFO] Starting Canopy-App simulation
2022-06-30 12:00:05 [INFO] Configuration loaded successfully
2022-06-30 12:00:10 [INFO] Input files validated
2022-06-30 12:00:15 [INFO] Beginning time integration
2022-06-30 12:01:00 [INFO] Time step 1/24 completed
2022-06-30 12:01:45 [INFO] Time step 2/24 completed
...
2022-06-30 12:02:30 [INFO] Simulation completed successfully
```

## Reading Output Files

### Python

Using netCDF4-python:

```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Open output file
ds = nc.Dataset('canopy_output.nc', 'r')

# Read variables
time = ds.variables['time'][:]
temp = ds.variables['temp'][:]
lai = ds.variables['lai'][:]

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(time, temp[:, 0, 0, 0])
plt.xlabel('Time (hours)')
plt.ylabel('Temperature (K)')
plt.title('Surface Temperature')
plt.show()

ds.close()
```

### NCO (NetCDF Operators)

Command-line tools for netCDF manipulation:

```bash
# Extract variable
ncks -v temp canopy_output.nc temp_only.nc

# Time average
ncwa -a time canopy_output.nc canopy_avg.nc

# Spatial subset
ncks -d lat,40.0,41.0 -d lon,-81.0,-80.0 canopy_output.nc subset.nc
```

### CDO (Climate Data Operators)

Advanced climate data operations:

```bash
# Time series statistics
cdo timstd canopy_output.nc temp_std.nc
cdo timmean canopy_output.nc temp_mean.nc

# Vertical interpolation
cdo intlevel,850,500,300 canopy_output.nc levels.nc
```

## Output Customization

### Namelist Options

Control output through namelist parameters:

```fortran
&CANOPY_OPTIONS
 ! Output format (1=netCDF, 2=text, 3=both)
 outfmt_opt = 1

 ! Output frequency (hours)
 output_freq = 1.0

 ! Variables to output
 output_vars = 'temp', 'qv', 'u', 'v', 'tke'

 ! Enable chemistry output
 output_chem = .true.

 ! Enable diagnostic output
 output_diag = .true.
/
```

### Variable Selection

Choose specific variables for output:

```fortran
&OUTPUT_CONTROL
 ! Meteorological variables
 out_temp = .true.
 out_qv = .true.
 out_wind = .true.
 out_tke = .true.

 ! Canopy variables
 out_ppfd = .true.
 out_tleaf = .true.
 out_gsw = .true.

 ! Chemical variables
 out_chem = .false.
/
```

## File Management

### Compression

Enable compression for large files:

```fortran
&CANOPY_OPTIONS
 compress_output = .true.
 compression_level = 6
/
```

### Chunking

Optimize file access patterns:

```fortran
&CANOPY_OPTIONS
 chunk_time = 24
 chunk_lev = 10
 chunk_lat = 32
 chunk_lon = 32
/
```

## Quality Control

The model performs output validation:

1. **Range checks**: Ensures values are physically reasonable
2. **Conservation checks**: Verifies mass and energy conservation
3. **Metadata validation**: Confirms proper units and attributes
4. **Format compliance**: Ensures CF convention compliance

For troubleshooting output issues, see the [Troubleshooting Guide](troubleshooting.md).
