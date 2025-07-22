# Examples & Tutorials

This section provides step-by-step examples and tutorials for common Canopy-App modeling scenarios.

## Basic Examples

### Example 1: Basic Forest Simulation

#### Scenario
Simulate a deciduous forest canopy for a summer day with meteorological forcing from GFS data.

#### Input Files Required

1. **Namelist file** (`namelist.canopy`)
2. **Meteorological data** (NetCDF format)
3. **Optional point data** for validation

#### Step 1: Prepare Input Data

```bash
# Create working directory
mkdir forest_example
cd forest_example

# Copy input files
cp ../input/gfs.t12z.20220701.sfcf000.canopy.nc .
cp ../input/namelist.canopy .
```

#### Step 2: Configure Namelist

Edit `namelist.canopy`:

```fortran
&filenames
    file_vars = 'gfs.t12z.20220701.sfcf000.canopy.nc'
    file_out = 'forest_output'
/

&userdefs
    ! Input format: 0=2D NetCDF, 1=1D text
    infmt_opt = 0

    ! Grid dimensions
    nlat = 361
    nlon = 720

    ! Time settings
    time_start = '2022-07-01-00:00:00.0000'
    ntime = 1
    time_intvl = 3600

    ! Canopy model vertical resolution
    modlays = 100
    modres = 0.5

    ! Enable physics components
    ifcanwind = .true.
    ifcanbio = .true.
    ifcanddepgas = .true.
    ifcanphot = .true.
    ifcaneddy = .true.

    ! Basic biogenic emission settings
    bio_cce = 0.21
    hist_opt = 0
/
```

#### Step 3: Run the Model

```bash
# Run Canopy-App
../canopy

# Check for successful completion
echo "Exit code: $?"
```

#### Step 4: Examine Output

```bash
# List output files
ls -la forest_output*

# If NetCDF tools are available
ncdump -h forest_output.nc | head -20
```

### Example 2: Point-Scale Simulation

#### Scenario
Run a single-point simulation for a specific location with text input files.

#### Step 1: Prepare Point Data

```bash
# Create point simulation directory
mkdir point_example
cd point_example

# Copy point input files
cp ../input/point_file_20220701.sfcf000.txt .
cp ../input/namelist.canopy .
```

#### Step 2: Configure for Point Mode

Edit `namelist.canopy`:

```fortran
&filenames
    file_vars = 'point_file_20220701.sfcf000.txt'
    file_out = 'point_output'
/

&userdefs
    ! Input format: 1=1D text
    infmt_opt = 1

    ! Single point
    nlat = 1
    nlon = 1

    ! Time settings
    time_start = '2022-07-01-00:00:00.0000'
    ntime = 1
    time_intvl = 3600

    ! Canopy model settings
    modlays = 50
    modres = 1.0

    ! Enable selected components
    ifcanwind = .true.
    ifcanbio = .true.
/
```

#### Step 3: Run and Analyze

```bash
# Run simulation
../canopy

# View text output
cat point_output.txt | head -20
```

## Advanced Examples

### Example 3: Multi-Day Biogenic Emissions with Historical Averaging

#### Scenario
Run a 10-day simulation to capture realistic biogenic emission patterns with historical temperature and light averaging.

#### Configuration

```fortran
&filenames
    file_vars = 'multi_day_input.nc'
    file_out = 'biogenic_10day'
/

&userdefs
    infmt_opt = 0
    nlat = 181
    nlon = 360

    ! 10-day simulation
    time_start = '2022-07-01-00:00:00.0000'
    ntime = 240  ! 10 days * 24 hours
    time_intvl = 3600  ! hourly

    modlays = 100
    modres = 0.5

    ! Biogenic emissions with advanced options
    ifcanbio = .true.
    bio_cce = 0.21
    hist_opt = 1  ! Enable historical averaging

    ! Soil moisture effects
    soim_opt = 0  ! Enable soil moisture response
    soild1 = 5.0
    soild2 = 25.0
    soild3 = 70.0
    soild4 = 150.0

    ! CO2 inhibition
    co2_opt = 0
    co2_set = 410.0  ! Current atmospheric CO2

    ! Leaf age effects
    leafage_opt = 0
    lai_tstep = 86400  ! Daily LAI updates
/
```

### Example 4: Fire Weather and Wind Adjustment Factor

#### Scenario
Calculate wind adjustment factors for wildfire applications during active fire conditions.

#### Configuration

```fortran
&filenames
    file_vars = 'fire_weather_input.nc'
    file_out = 'fire_waf_output'
/

&userdefs
    infmt_opt = 0
    nlat = 100
    nlon = 100

    time_start = '2022-08-15-12:00:00.0000'
    ntime = 24  ! 24 hours
    time_intvl = 3600

    ! Fine resolution for WAF calculations
    modlays = 200
    modres = 0.25  ! 25 cm resolution

    ! Wind and fire components
    ifcanwind = .true.
    ifcanwaf = .true.

    ! Wind settings
    href_opt = 0
    href_set = 10.0  ! 10m reference height
    z0ghc = 0.1
    rsl_opt = 0
    lambdars = 2.0

    ! Fire/WAF settings
    dx_opt = 1
    dx_set = 1000.0  ! 1km grid resolution
    flameh_opt = 0  ! Calculate from FRP
    flameh_cal = 1  ! Crown scorch method
    frp_fac = 1.0
/
```

### Example 5: Urban Dry Deposition

#### Scenario
Model gas dry deposition in an urban environment with building surfaces.

#### Configuration

```fortran
&filenames
    file_vars = 'urban_input.nc'
    file_out = 'urban_drydep'
/

&userdefs
    infmt_opt = 0
    nlat = 50
    nlon = 50

    time_start = '2022-06-15-06:00:00.0000'
    ntime = 48  ! 2 days
    time_intvl = 3600

    modlays = 80
    modres = 0.5

    ! Dry deposition
    ifcanddepgas = .true.

    ! Gas chemistry settings
    chemmechgas_opt = 0  ! RACM2
    chemmechgas_tot = 31
    ddepspecgas_opt = 0  ! All species

    ! Urban-specific settings
    gamma_set = 5.0e-5  ! Building surface reactivity
    Ramin_set = 10.0    ! Minimum aerodynamic resistance
    hyblev1 = 20.0      ! First model level height

    ! Snow/ice thresholds
    snowc_set = 50.0
    icec_set = 50.0
/
```

## Analysis and Visualization Examples

### Python Analysis Script

```python
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Load Canopy-App output
ds = xr.open_dataset('forest_output.nc')

# Plot vertical wind profile
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(ds.ws[0, 0, 0, :], ds.z, 'b-', linewidth=2)
ax.set_xlabel('Wind Speed (m/s)')
ax.set_ylabel('Height (m)')
ax.set_title('Canopy Wind Profile')
ax.grid(True)
plt.show()

# Plot biogenic emissions profile
if 'emi_isop' in ds.variables:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(ds.emi_isop[0, 0, 0, :] * 1e9, ds.z, 'g-', linewidth=2)
    ax.set_xlabel('Isoprene Emission (ng/m³/s)')
    ax.set_ylabel('Height (m)')
    ax.set_title('Vertical Isoprene Emission Profile')
    ax.grid(True)
    plt.show()
```

### Validation Against Observations

```python
# Load observation data
obs_data = np.loadtxt('tower_observations.txt')
obs_height = obs_data[:, 0]
obs_wind = obs_data[:, 1]

# Load model output
model_wind = ds.ws[0, 0, 0, :].values
model_height = ds.z.values

# Interpolate model to observation heights
from scipy.interpolate import interp1d
f = interp1d(model_height, model_wind, kind='linear')
model_interp = f(obs_height)

# Calculate statistics
rmse = np.sqrt(np.mean((model_interp - obs_wind)**2))
corr = np.corrcoef(model_interp, obs_wind)[0, 1]

print(f"RMSE: {rmse:.3f} m/s")
print(f"Correlation: {corr:.3f}")
```

## Common Issues and Solutions

### Issue 1: Model Fails to Start

**Problem**: Error message "NetCDF file not found"

**Solution**:
```bash
# Check file existence and permissions
ls -la input_file.nc
# Verify NetCDF file integrity
ncdump -h input_file.nc
```

### Issue 2: Unrealistic Wind Speeds

**Problem**: Wind speeds are too high or too low

**Solutions**:
- Check `href_set` value (recommend 10m)
- Verify `z0ghc` ratio (typical: 0.1)
- Ensure `modres` is appropriate for canopy height

### Issue 3: No Biogenic Emissions

**Problem**: All emission values are zero

**Solutions**:
- Verify LAI > 0 in input data
- Check temperature is reasonable (> 0°C)
- Ensure `bio_cce` is set (recommend 0.21)
- Check if `co2_opt` is inhibiting emissions

### Issue 4: WAF Calculation Errors

**Problem**: Wind adjustment factor values are unrealistic

**Solutions**:
- Use fine vertical resolution (`modres <= 0.5m`)
- Check flame height calculations
- Verify FRP input data quality
- Ensure grid resolution `dx_set` is appropriate

## Performance Tips

### For Large Domains
- Use coarser vertical resolution (`modres = 1.0-2.0m`)
- Reduce number of vertical layers (`modlays = 50-75`)
- Disable unnecessary components

### For Long Simulations
- Enable historical averaging (`hist_opt = 1`)
- Use appropriate spin-up period (240 hours recommended)
- Monitor memory usage for large grids

### For High-Resolution Studies
- Use fine vertical resolution (`modres = 0.25-0.5m`)
- Increase vertical layers (`modlays = 100-200`)
- Ensure adequate computational resources

## Next Steps

- **[Namelist Reference](namelist-reference.md)** - Complete parameter documentation
- **[Troubleshooting](troubleshooting.md)** - Detailed problem-solving guide
- **[Science Documentation](../science/model-description.md)** - Understanding the physics
- **[API Reference](../api/overview.md)** - Code documentation
