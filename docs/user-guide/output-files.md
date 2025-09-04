# Output Files

The Canopy-App model generates various output files containing simulation results. This section describes the format and content of each output file type.

## Output File Types

### NetCDF Output Files

The model generates netCDF files with comprehensive simulation results:

#### Primary Output File: `canopy_output.nc`

Contains the main simulation variables:

**Dimensions:**
- `time`: Number of time steps
- `lev`: Number of vertical levels
- `lat`: Number of latitude points
- `lon`: Number of longitude points

**Example Canopy-Process Variables:**

| Variable | Description | Units | Dimensions |
|-----------|---------------------------------------|-------------------|--------------------|
| `time`    | Time coordinate                       | hours since start | (time)             |
| `lev`     | Vertical level coordinate             | m                 | (lev)              |
| `lat`     | Latitude coordinate                   | degrees_north     | (lat)              |
| `lon`     | Longitude coordinate                  | degrees_east      | (lon)              |
| `canwind` | Canopy winds                          | m/s               | (time,lev,lat,lon) |
| `kz`      | Canopy eddy diffusivities             | kg/kg             | (time,lev,lat,lon) |
| `rjcf`    | Canopy photolysis attenuation factors |                   | (time,lev,lat,lon) |
| `emi_isop`| Biogenic emissions                    | kg/m3 s           | (time,lev,lat,lon) |
| `ddep_o3` | Ozone dry deposition                  | cm/s              | (time,lev,lat,lon) |
|---------------------------------------------------------------------------------------------

#### Example Canopy-Specific Variables

| Variable | Description | Units | Dimensions |
| `canheight` | Canopy height | m | (time,lat,lon) |
| `lad` | Leaf area density | m²/m3 | (time,lev,lat,lon) |
| `z0_h` | Ratio of surface roughness length to canheight |  | (time,lev,lat,lon) |
| `gsw_shade` | Shaded stomatal conductance | mol/m²/s | (time,lev,lat,lon) |

### Text Output Files

For simple analysis, the model can output text files:

#### Point Output: `point_file_YYYYMMDD.sfcfXXX.txt`

Contains time series data for single points:

```
# Example Canopy-App Point Output for Canopy Winds
#    time stamp: 2022-07-01-11:00:00.0000
#    reference height, h:   10.0 m
#    number of model layers:    100
#    lat      lon  height (m)  LAD (m2 m-3)       ws (m s-1)
      34.97   270.00      0.00        0.00  0.0000000E+00
      34.97   270.00      0.50        0.08  5.9679881E-02
      34.97   270.00      1.00        0.13  7.8142006E-02
    ...
```

#### Example 3D Canopy Variables for Input PAVD: `point_file_canvars_YYYYMMDD.sfcfXXX.txt`

<!-- ```
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

<!-- ## Chemical Output (if enabled)

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
| `so2` | SO₂ concentration | ppb | (time,lev,lat,lon) | -->

### Emission Rates

Example biogenic emission rates:

| Variable | Description | Units | Dimensions |
|----------|-------------|-------|------------|
| `emi_isop` | Isoprene emission rate | kg/m²/s | (time,lev,lat,lon) |
| `emi_mono` | Monoterpene emission rate | kg/m²/s | (time,lev,lat,lon) |
| `emi_sesq` | Sesquiterpene emission rate | kg/m²/s | (time,lev,lat,lon) |
| `emi_ovoc` | Other VOC emission rate | kg/m²/s | (time,lev,lat,lon) |

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
 ! Output format (1=netCDF, 2=text)
 outfmt_opt = 1

/
```


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

## Quality Control

The model performs output validation:

1. **Range checks**: Ensures values are physically reasonable
2. **Conservation checks**: Verifies mass and energy conservation
3. **Metadata validation**: Confirms proper units and attributes
4. **Format compliance**: Ensures CF convention compliance

For troubleshooting output issues, see the [Troubleshooting Guide](troubleshooting.md).
