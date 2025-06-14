# Configuration Guide

Learn how to configure Canopy-App for your specific modeling needs through the Fortran namelist system.

## Overview

Canopy-App uses a Fortran namelist file for all configuration settings. The main configuration file is `input/namelist.canopy` and contains several namelist groups that control input/output, model physics, and component-specific options.

## Namelist Structure

The namelist file is organized into logical groups:

```fortran
&filenames
    ! Input/Output file specifications
/

&userdefs
    ! Core model configuration and physics options
/
```

## Example Configuration

Here's a basic example of `input/namelist.canopy`:

```fortran
&filenames
    file_vars = 'input/gfs.t12z.20220701.sfcf000.canopy.nc'
    file_out = 'output_test'
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
    
    ! Physics components
    ifcanwind = .true.
    ifcanbio = .true.
    ifcanddepgas = .true.
    ifcanphot = .true.
    ifcaneddy = .true.
/
```

## Complete Namelist Options

### Input/Output Configuration

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `infmt_opt` | integer | Input format: `0`=2D NetCDF, `1`=1D text | `0` |
| `nlat` | integer | Number of latitude cells | - |
| `nlon` | integer | Number of longitude cells | - |
| `file_vars` | string | Input file path | - |
| `file_out` | string | Output file prefix | - |

### Time Configuration

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `time_start` | string | Start time (YYYY-MM-DD-HH:MM:SS.SSSS) | - |
| `time_end` | string | End time (YYYY-MM-DD-HH:MM:SS.SSSS) | - |
| `ntime` | integer | Number of time steps | - |
| `time_intvl` | integer | Time interval in seconds | - |

### Vegetation and Land Use

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `lu_opt` | integer | Land use type: `0`=VIIRS 17-cat, `1`=MODIS-IGBP 20-cat | `0` |
| `var3d_opt` | integer | Use 3D variables: `0`=off, `1`=on | `0` |
| `var3d_set` | integer | Number of 3D input levels | `14` |
| `pavd_opt` | integer | Use GEDI 3D PAVD profiles: `0`=off, `1`=on | `0` |
| `pavd_set` | real | Latitude threshold for PAVD profiles | `52.0` |

### Canopy Model Structure

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `modlays` | integer | Number of model layers | - |
| `modres` | real | Vertical resolution (m) | - |
| `lai_thresh` | real | LAI threshold for contiguous canopy (m²/m²) | - |
| `cf_thresh` | real | Canopy fraction threshold | - |
| `ch_thresh` | real | Canopy height threshold (m) | - |

### Vegetation Type Settings

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `can_opt` | integer | Use namelist values: `0`=input data, `1`=namelist | `0` |
| `can_chset` | real | Constant canopy height (m) | - |
| `can_cfset` | real | Constant canopy fraction | - |
| `can_laiset` | real | Constant canopy LAI | - |
| `ssg_opt` | integer | Shrub/savanna/grass option | `0` |
| `ssg_chset` | real | SSG height (m) | - |
| `ssg_cfset` | real | SSG fraction | - |
| `ssg_laiset` | real | SSG LAI | - |
| `crop_opt` | integer | Crop vegetation option | `0` |
| `crop_chset` | real | Crop height (m) | - |
| `crop_cfset` | real | Crop fraction | - |
| `crop_laiset` | real | Crop LAI | - |

## Physics Component Configuration

### Wind Calculations

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanwind` | logical | Enable canopy wind calculations | `.FALSE.` |
| `href_opt` | integer | Reference height: `0`=namelist, `1`=file | `0` |
| `href_set` | real | Reference height above canopy (m) | - |
| `z0ghc` | real | Ground roughness to canopy height ratio | - |
| `rsl_opt` | integer | Roughness sublayer: `0`=MOST, `1`=unified RSL | `0` |
| `lambdars` | real | RSL influence factor | - |
| `pai_opt` | integer | PAI calculation method (0-3) | `0` |
| `pai_set` | real | User-set PAI value | `4.0` |
| `z0_opt` | integer | Surface roughness: `0`=input, `1`=vegetation-dependent | `0` |

### Wind Adjustment Factor (WAF)

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanwaf` | logical | Enable WAF calculations | `.FALSE.` |
| `dx_opt` | integer | Grid resolution: `0`=calculate, `1`=user-set | `0` |
| `dx_set` | real | Grid resolution (m) | - |
| `flameh_opt` | integer | Flame height method (0-5) | - |
| `flameh_cal` | integer | FRP calculation method (0-1) | - |
| `flameh_set` | real | Flame height (m) or fraction | - |
| `frp_fac` | real | FRP tuning factor | `1.0` |

### Vertical Diffusion

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcaneddy` | logical | Enable eddy diffusivity calculations | `.FALSE.` |

### Photolysis

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanphot` | logical | Enable photolysis attenuation | `.FALSE.` |

### Biogenic Emissions

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanbio` | logical | Enable biogenic emissions | `.FALSE.` |
| `bio_cce` | real | Canopy environment coefficient | `0.21` |
| `biospec_opt` | integer | Species output: `0`=all, `1-19`=single species | `0` |
| `biovert_opt` | integer | Vertical summing option (0-3) | `0` |
| `loss_opt` | integer | Canopy loss factor (0-2) | `0` |
| `loss_set` | real | Constant loss factor | `0.96` |
| `loss_ind` | integer | Species for loss factor | `0` |
| `lifetime` | real | Chemical lifetime (s) | `3600` |
| `co2_opt` | integer | CO₂ inhibition method (0-2) | `0` |
| `co2_set` | real | CO₂ concentration (ppmv) | - |
| `leafage_opt` | integer | Leaf age response: `0`=on, `1`=off | `1` |
| `lai_tstep` | integer | LAI data interval (s) | - |
| `hist_opt` | integer | Historical averaging: `0`=off, `1`=on | `0` |
| `soim_opt` | integer | Soil moisture response: `0`=on, `1`=off | `1` |
| `soild1-4` | real | Soil layer depths (cm) | `5.0, 25.0, 70.0, 150.0` |
| `aq_opt` | integer | Air quality stress (0-2) | `2` |
| `w126_set` | real | Ozone W126 value (ppm-hours) | - |
| `ht_opt` | integer | High temperature stress: `0`=on, `1`=off | `1` |
| `lt_opt` | integer | Low temperature stress: `0`=on, `1`=off | `1` |
| `hw_opt` | integer | High wind stress: `0`=on, `1`=off | `1` |

### Gas Dry Deposition

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanddepgas` | logical | Enable gas dry deposition | `.FALSE.` |
| `ddepspecgas_opt` | integer | Species output: `0`=all, `1-31`=single | `0` |
| `chemmechgas_opt` | integer | Chemical mechanism: `0`=RACM2 | `0` |
| `chemmechgas_tot` | integer | Total species in mechanism | `31` |
| `hyblev1` | real | Height of first hybrid level (m) | `20.0` |
| `snowc_set` | real | Snow cover threshold (%) | `50.0` |
| `icec_set` | real | Ice cover threshold (%) | `50.0` |
| `gamma_set` | real | Building surface reaction probability | `5.0e-5` |
| `Ramin_set` | real | Minimum aerodynamic resistance (s/m) | `10.0` |

## Configuration Tips

### For Beginners
1. **Start Simple**: Enable one physics component at a time
2. **Use Defaults**: Most parameters have sensible defaults
3. **Check Units**: Pay attention to units in variable descriptions
4. **Validate Settings**: Use reasonable values for your domain

### For Advanced Users
1. **Historical Averaging**: Use `hist_opt=1` for biogenic emissions with >24hr runs
2. **Spin-up**: Run 10 days (240 hrs) for optimal biogenic emission analysis
3. **Resolution**: Match `modres` and `modlays` to your canopy heights
4. **WAF Calculations**: Use fine resolution (`modres <= 0.5 m`) for WAF

### Common Configurations

#### Basic Wind and Emissions
```fortran
&userdefs
    ifcanwind = .true.
    ifcanbio = .true.
    bio_cce = 0.21
    hist_opt = 0  ! Use 1 for runs >24 hours
/
```

#### Full Physics Suite
```fortran
&userdefs
    ifcanwind = .true.
    ifcanwaf = .true.
    ifcaneddy = .true.
    ifcanphot = .true.
    ifcanbio = .true.
    ifcanddepgas = .true.
    
    ! Biogenic emissions with historical averaging
    hist_opt = 1
    soim_opt = 0  ! Enable soil moisture effects
/
```

## Validation and Testing

After configuring your namelist:

1. **Check Syntax**: Ensure proper Fortran namelist format
2. **Validate Ranges**: Check that values are within reasonable bounds
3. **Test Run**: Start with a short simulation to verify settings
4. **Review Output**: Examine output variables for expected ranges

## Next Steps

- **[Run the Model](quickstart.md)**: Execute your configured simulation
- **[Examples](../examples/basic.md)**: See practical configuration examples
- **[Science Guide](../science/model-description.md)**: Understand the physics behind the options
|--------|-------|-------------|
| `infmt_opt` | 1 | NetCDF input format |
| `infmt_opt` | 2 | Text input format |
| `outfmt_opt` | 1 | NetCDF output format |
| `outfmt_opt` | 2 | Text output format |

## Model Options

### Physics Switches

```fortran
&canopy_options
    ! Core physics
    opt_canmet  = 1             ! Canopy meteorology (0=off, 1=on)
    opt_bioem   = 1             ! Biogenic emissions (0=off, 1=on)
    opt_drydep  = 1             ! Dry deposition (0=off, 1=on)
    opt_phot    = 1             ! Photolysis rates (0=off, 1=on)

    ! Radiation options
    opt_rad     = 1             ! Radiation transfer (0=off, 1=on)
    opt_solarzen = 1            ! Solar zenith angle calc (0=off, 1=on)

    ! Chemical options
    opt_chem    = 0             ! Chemistry (0=off, 1=on)
/
```

### Advanced Options

```fortran
&canopy_physics
    ! Canopy structure
    z_canbot    = 0.0           ! Canopy bottom height (m)
    z_cantop    = 20.0          ! Canopy top height (m)
    ncanlevs    = 10            ! Number of canopy levels

    ! Meteorology
    stability_opt = 1           ! Stability correction (1=on, 0=off)

    ! Emissions
    bioem_opt   = 1             ! Biogenic emission algorithm
    temp_opt    = 1             ! Temperature dependence option

    ! Deposition
    drydep_opt  = 1             ! Dry deposition algorithm
    stom_opt    = 1             ! Stomatal resistance option
/
```

## Configuration Examples

### Basic Research Setup

```fortran
&canopy_inputs
    in_date = '20220701'
    in_time = '12'
    file_in = 'gfs_input.nc'
    file_out = 'canopy_research.nc'
    infmt_opt = 1
    outfmt_opt = 1
/

&canopy_options
    opt_canmet = 1
    opt_bioem = 1
    opt_drydep = 1
    opt_phot = 1
    opt_rad = 1
/
```

### Production Run Setup

```fortran
&canopy_inputs
    in_date = '20220701'
    in_time = '00'
    file_in = 'operational_input.nc'
    file_out = 'operational_output.nc'
    file_pnt = 'point_output.txt'
    infmt_opt = 1
    outfmt_opt = 1
/

&canopy_options
    opt_canmet = 1
    opt_bioem = 1
    opt_drydep = 1
    opt_phot = 1
    opt_rad = 1
    opt_solarzen = 1
/

&canopy_physics
    z_cantop = 25.0
    ncanlevs = 15
    stability_opt = 1
    bioem_opt = 2
    drydep_opt = 2
/
```

### Sensitivity Study Setup

```fortran
&canopy_inputs
    in_date = '20220701'
    in_time = '12'
    file_in = 'sensitivity_input.nc'
    file_out = 'sensitivity_test.nc'
/

&canopy_options
    opt_canmet = 1
    opt_bioem = 1
    opt_drydep = 0    ! Turn off dry deposition
    opt_phot = 1
    opt_rad = 1
/

&canopy_physics
    z_cantop = 15.0   ! Lower canopy height
    ncanlevs = 8      ! Fewer levels for speed
    bioem_opt = 1     ! Basic emission algorithm
/
```

## Validation and Testing

### Configuration Validation

The model performs several validation checks:

- Date and time format validation
- File existence checks
- Parameter range validation
- Physics option compatibility

### Testing Your Configuration

```bash
# Test configuration without full run
./canopy_app --check-config

# Run with verbose output
./canopy_app --verbose

# Dry run (parse inputs only)
./canopy_app --dry-run
```

## Common Configuration Issues

!!! error "File Not Found"
    **Problem**: Input files not found
    **Solution**: Check file paths are relative to run directory

!!! warning "Parameter Out of Range"
    **Problem**: Physics parameters outside valid ranges
    **Solution**: Check parameter documentation and valid ranges

!!! tip "Performance Tuning"
    **Tip**: Adjust `ncanlevs` based on computational resources and accuracy needs

## Advanced Configuration

### Environment Variables

```bash
# Set NetCDF library path
export NETCDF=/usr/local/netcdf

# Set number of OpenMP threads
export OMP_NUM_THREADS=4

# Set stack size for large arrays
ulimit -s unlimited
```

### Custom Physics Parameters

For advanced users, physics parameters can be modified in the source code modules:

- `canopy_const_mod.F90` - Physical constants
- `canopy_bioparm_mod.F90` - Biogenic emission parameters
- `canopy_canopts_mod.F90` - Model options and defaults

## Next Steps

- [Learn about input file formats](../user-guide/input-files.md)
- [Understand model physics](../science/model-description.md)
- [Explore examples](../examples/basic.md)
- [Check troubleshooting guide](../user-guide/troubleshooting.md)
