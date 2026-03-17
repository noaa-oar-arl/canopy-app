# Namelist Reference

Complete reference guide for all Canopy-App namelist options and parameters.

## Overview

Canopy-App uses Fortran namelist files for configuration. The main configuration file is [`input/namelist.canopy`](https://github.com/noaa-oar-arl/canopy-app/blob/main/input/namelist.canopy) and contains two primary namelist groups: `&filenames` and `&userdefs`.

## Example Namelist File

```fortran
&filenames
    file_vars = 'input/gfs.t12z.20220701.sfcf000.canopy.nc'
    file_out = 'canopy_output'
/

&userdefs
    ! Input format options
    infmt_opt = 0
    nlat = 361
    nlon = 720

    ! Time configuration
    time_start = '2022-07-01-00:00:00.0000'
    time_end = '2022-07-01-23:00:00.0000'
    ntime = 24
    time_intvl = 3600

    ! Canopy model structure
    modlays = 100
    modres = 0.5

    ! Physics components
    ifcanwind = .true.
    ifcanbio = .true.
    ifcanddepgas = .true.
    ifcanphot = .true.
    ifcaneddy = .true.

    ! Biogenic emission settings
    bio_cce = 0.21
    hist_opt = 1
    soim_opt = 0
/
```

## Complete Namelist Options Reference

### Input Model Format Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `infmt_opt` | integer | Input format: `0`=2D NetCDF, `1`=1D text | `0` |

### Input Model Grid Sizes

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `nlat` | integer | Number of latitude cells (must match input file) | - |
| `nlon` | integer | Number of longitude cells (must match input file) | - |

### Input Model Run Times and Interval

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `time_start` | string | Start time (YYYY-MM-DD-HH:MM:SS.SSSS) | - |
| `time_end` | string | End time (YYYY-MM-DD-HH:MM:SS.SSSS) | - |
| `ntime` | integer | Number of time steps | - |
| `time_intvl` | integer | Time interval in seconds (e.g., 3600 for hourly) | - |

### Canopy Model Vegetation/Land Use Input Dataset Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `lu_opt` | integer | Land use type: `0`=VIIRS 17-cat, `1`=MODIS-IGBP 20-cat | `0` |

### Input Model 3D NetCDF Variable Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `var3d_opt` | integer | Use 3D variables: `0`=off, `1`=on | `0` |
| `var3d_set` | integer | Number of 3D input levels (max 14 for text files) | `14` |

### Options to Use Observed PAVD Profiles and Latitude Threshold

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `pavd_opt` | integer | Use GEDI 3D PAVD profiles: `0`=off, `1`=on | `0` |
| `pavd_set` | real | Latitude threshold for PAVD profiles (degrees) | `52.0` |

### Canopy Model Vertical Layers

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `modlays` | integer | Number of model layers (above and below canopy) | - |
| `modres` | real | Vertical resolution (m) | - |

!!! tip "Layer Configuration"
    Strongly recommend adjusting `modlays` and `modres` to maintain canopy model column extension above tallest canopies. For a 50m column simulation, you could use:
    - 1000 layers @ 0.05m resolution
    - 100 layers @ 0.5m resolution
    - 50 layers @ 1.0m resolution

### Contiguous Canopy Model Thresholds

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `lai_thresh` | real | LAI threshold for contiguous canopy (m²/m²) | - |
| `cf_thresh` | real | Canopy fraction threshold for contiguous canopy | - |
| `ch_thresh` | real | Canopy height threshold for contiguous canopy (m) | - |

!!! note
    These thresholds only apply for valid vegetated land use types (forests, shrub/savanna/grass, crops, wetlands).

### Canopy Crop and Shrub/Savanna/Grass Extension Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `can_opt` | integer | Use input data (`0`) or user-set values (`1`) | `0` |
| `can_chset` | real | Constant canopy height (m), used if `can_opt=1` | - |
| `can_cfset` | real | Constant canopy fraction, used if `can_opt=1` | - |
| `can_laiset` | real | Constant canopy LAI, used if `can_opt=1` | - |
| `ssg_opt` | integer | Shrub/savanna/grass option: input data (`0`) or user-set (`1`) | `0` |
| `ssg_chset` | real | Constant SSG height (m), used if `ssg_opt=1` | `0.5-1.0` |
| `ssg_cfset` | real | Constant SSG fraction, used if `ssg_opt=1` | - |
| `ssg_laiset` | real | Constant SSG LAI, used if `ssg_opt=1` | - |
| `crop_opt` | integer | Crop option: input data (`0`) or user-set (`1`) | `0` |
| `crop_chset` | real | Constant crop height (m), used if `crop_opt=1` | - |
| `crop_cfset` | real | Constant crop fraction, used if `crop_opt=1` | - |
| `crop_laiset` | real | Constant crop LAI, used if `crop_opt=1` | - |

### Canopy Physics and Wind-Specific Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanwind` | logical | Enable canopy wind calculations | `.FALSE.` |
| `href_opt` | integer | Use namelist (`0`) or file array (`1`) for reference height | `0` |
| `href_set` | real | Reference height above canopy (m) | `10.0` |
| `z0ghc` | real | Ratio of ground roughness to canopy height | - |
| `rsl_opt` | integer | RSL effects: MOST (`0`) or unified RSL (`1`) | `0` |
| `lambdars` | real | RSL influence factor (with `rsl_opt=0`) | - |
| `pai_opt` | integer | PAI calculation method (`0-3`) | `0` |
| `pai_set` | real | User-set PAI value (used if `pai_opt=3`) | `4.0` |
| `z0_opt` | integer | Use model input (`0`) or vegetation-dependent z0 (`1`) | `0` |

!!! warning
    If `href_set` becomes small and approaches z0, only sub-canopy wind profile is calculated. Recommend `href_set = 10m`.

### Canopy Fire/WAF-Specific Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanwaf` | logical | Enable canopy WAF calculations | `.FALSE.` |
| `dx_opt` | integer | Grid resolution: calculated (`0`) or user-set (`1`) | `0` |
| `dx_set` | real | Grid resolution (m), used if `dx_opt=1` | - |
| `flameh_opt` | integer | Flame height calculation method (`0-5`) | `0` |
| `flameh_cal` | integer | FRP-based flame height method (`0-1`) | `0` |
| `flameh_set` | real | User-set flame height (m) or fraction of canopy height | - |
| `frp_fac` | real | Tuning factor for FRP in flame height calculation | `1.0` |

!!! warning
    If `modres >> flameh`, errors in WAF calculation will occur. Use fine `modres` (≤0.5m) compared to average flame heights (~1.0m).

### Canopy Eddy Diffusivity-Specific Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcaneddy` | logical | Enable canopy eddy Kz calculations | `.FALSE.` |

### Canopy Radiation/Photolysis-Specific Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanphot` | logical | Enable canopy photolysis calculations | `.FALSE.` |

### Canopy Biogenic Emissions-Specific Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanbio` | logical | Enable canopy biogenic emissions | `.FALSE.` |
| `bio_cce` | real | MEGAN canopy environment coefficient | `0.21` |
| `biospec_opt` | integer | Species output: all (`0`) or single species (`1-19`) | `0` |
| `biovert_opt` | integer | Vertical summing option (`0-3`) | `0` |
| `loss_opt` | integer | Canopy loss factor option (`0-2`) | `0` |
| `loss_set` | real | Constant loss factor (used if `loss_opt=2`) | `0.96` |
| `loss_ind` | integer | Apply loss to all species (`0`) or specific species | `0` |
| `lifetime` | real | Chemical lifetime in seconds (used if `loss_opt=1`) | `3600` |
| `co2_opt` | integer | CO₂ inhibition method: Possell & Hewitt (`0`), Wilkinson (`1`), off (`2`) | `0` |
| `co2_set` | real | Atmospheric CO₂ concentration (ppmv) | - |

For complete details on all biogenic emission options including leaf age, historical averaging, soil moisture, and environmental stress factors, see the full reference in the [README.md](https://github.com/noaa-oar-arl/canopy-app/blob/main/README.md#table-3-current-user-namelist-options).

### Canopy Gas Dry Deposition-Specific Options

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanddepgas` | logical | Enable canopy gas dry deposition | `.FALSE.` |
| `ddepspecgas_opt` | integer | Species output: all (`0`) or single species (`1-31`) | `0` |
| `chemmechgas_opt` | integer | Chemical mechanism: RACM2 (`0`) | `0` |
| `chemmechgas_tot` | integer | Total species in mechanism | `31` |
| `hyblev1` | real | Height of 1st hybrid model layer (m) | `20.0` |
| `snowc_set` | real | Snow cover threshold (%) | `50.0` |
| `icec_set` | real | Ice cover threshold (%) | `50.0` |
| `gamma_set` | real | Building surface reaction probability | `5.0e-5` |
| `Ramin_set` | real | Minimum aerodynamic resistance (s/m) | `10.0` |

### Soil NO Emissions and Canopy Reduction Factor (CRF)

| Option      | Type    | Description                                              | Default |
|-------------|---------|----------------------------------------------------------|---------|
| `soilno_opt`| integer | Soil NO emissions: `0`=off, `1`=use BDSNP model         | `1`     |
| `crf_opt`   | integer | Canopy reduction factor: `0`=default bulk Beer's law, `1`=LAD-resolved Beer's law, `2`=BDSNP NO₂ deposition (LAI), `3`=BDSNP NO₂ deposition (LAD). Options 2 and 3 require `ifcanwind=.TRUE.` | `1`     |
| `crf_set`   | real    | User override for CRF (0=auto, else fixed value)        | `0.0`   |
| `fert_frac` | real    | Fraction of applied fertilizer N emitted as soil NO. `0.01`=1% (Steinkamp & Lawrence 2011), `0.025`=2.5% (Hudman et al. 2012) | `0.01` |
| `soiltemp_opt` | integer | Soil temperature response: `0`=Q10 unbounded, `1`=YL95/BDSNP (saturates at 30°C, zero below 0°C), `2`=Wang et al. (2021) cubic for T>20°C | `0` |

These options control the new soil NO emissions module and the method for calculating the canopy reduction factor (CRF). Set `soilno_opt=1` to enable the BDSNP-based soil NO emissions. CRF options: `crf_opt=0` uses bulk `exp(−0.5×LAI)`; `crf_opt=1` uses the vertically resolved LAD profile with Beer's law; `crf_opt=2` uses BDSNP physics-based NO₂ deposition velocity with bulk LAI; `crf_opt=3` uses BDSNP NO₂ deposition velocity with per-layer LAD (most physically complete). Options 2 and 3 require `ifcanwind=.TRUE.` to provide in-canopy meteorological profiles. If `crf_set` is nonzero, it overrides the computed CRF value. The `fert_frac` parameter controls the fraction of applied fertilizer nitrogen (from NPKGRIDS) that is emitted as soil NO. The `soiltemp_opt` parameter selects the soil temperature response function: `0` uses the standard Q10 exponential without bounds; `1` uses the YL95/BDSNP formulation where the response saturates at 30°C and drops to zero below 0°C; `2` uses the Wang et al. (2021) observation-based cubic function (−0.009T³+0.837T²−22.527T+196.149) above 20°C, which allows emissions to continue increasing in warm soils (20–40°C) rather than capping at 30°C.

## Configuration Tips

### Basic Setup
- Start with default values and enable only needed physics components
- Ensure grid dimensions (`nlat`, `nlon`) match your input files
- Set appropriate vertical resolution (`modres`) and layers (`modlays`)

### Physics Components
- Enable components based on your application:
  - `ifcanwind`: Wind profiles, wildfire applications
  - `ifcanbio`: Biogenic emissions
  - `ifcanddepgas`: Gas dry deposition
  - `ifcanphot`: Photolysis attenuation
  - `ifcaneddy`: Vertical diffusion

### Performance Considerations
- Higher vertical resolution (`modres`) provides more accuracy but increases computation time
- Use historical averaging (`hist_opt=1`) for longer simulations to improve biogenic emission calculations
- Consider enabling soil moisture response (`soim_opt=0`) for realistic biogenic emissions

### Common Issues
- Ensure `href_set` ≥ 10m to avoid calculation issues
- Match time settings with your input data frequency
- Verify land use type compatibility with your input data
- Check that flame height settings are appropriate for your grid resolution

## Example Configurations

### Basic Forest Simulation
```fortran
&userdefs
    ! Basic setup
    infmt_opt = 0
    nlat = 361, nlon = 720
    modlays = 100, modres = 0.5

    ! Enable core physics
    ifcanwind = .true.
    ifcanbio = .true.
    ifcanphot = .true.

    ! Basic biogenic settings
    bio_cce = 0.21
    hist_opt = 0
/
```

### Advanced Biogenic Emissions
```fortran
&userdefs
    ! Detailed biogenic setup
    ifcanbio = .true.
    bio_cce = 0.21
    hist_opt = 1
    soim_opt = 0

    ! Environmental factors
    aq_opt = 0
    co2_opt = 0
    co2_set = 415.0
/
```

### Wildfire Application
```fortran
&userdefs
    ! Wind and fire setup
    ifcanwind = .true.
    ifcanwaf = .true.

    ! Fine resolution for flame height
    modres = 0.25
    modlays = 200

    ! Fire-specific settings
    flameh_opt = 0
    flameh_cal = 1
    frp_fac = 1.0
/
```

### Canopy Model Vertical Structure

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `modlays` | integer | Number of model layers (above/below canopy) | - |
| `modres` | real | Vertical resolution (m) | - |
| `lai_thresh` | real | LAI threshold for contiguous canopy (m²/m²) | - |
| `cf_thresh` | real | Canopy fraction threshold | - |
| `ch_thresh` | real | Canopy height threshold (m) | - |

### Vegetation Type Override Options

#### General Canopy Settings
| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `can_opt` | integer | Use namelist values: `0`=input data, `1`=namelist values when missing | `0` |
| `can_chset` | real | Constant canopy height (m) | - |
| `can_cfset` | real | Constant canopy fraction | - |
| `can_laiset` | real | Constant canopy LAI | - |

#### Shrub/Savanna/Grass (SSG) Settings
| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ssg_opt` | integer | Use namelist SSG values when missing: `0`=off, `1`=on | `0` |
| `ssg_chset` | real | SSG height (m) - recommend 0.5-1.0m | - |
| `ssg_cfset` | real | SSG fraction | - |
| `ssg_laiset` | real | SSG LAI | - |

#### Crop Settings
| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `crop_opt` | integer | Use namelist crop values when missing: `0`=off, `1`=on | `0` |
| `crop_chset` | real | Crop height (m) | - |
| `crop_cfset` | real | Crop fraction | - |
| `crop_laiset` | real | Crop LAI | - |

## Physics Component Options

### Wind Calculations

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanwind` | logical | Enable canopy wind calculations | `.FALSE.` |
| `ifcanaeroddep` | logical | Enable sub-canopy aerosol dry deposition (Katul et al. 2010) | `.FALSE.` |
| `aeroddep_diam` | real | Aerosol particle diameter (m) | `1.0e-6` |
| `aeroddep_rho` | real | Aerosol particle density (kg/m^3) | `1500.0` |
| `href_opt` | integer | Reference height: `0`=use `href_set`, `1`=from file | `0` |
| `href_set` | real | Reference height above canopy (m) - recommend 10m | - |
| `z0ghc` | real | Ground roughness to canopy height ratio | - |
| `rsl_opt` | integer | Roughness sublayer: `0`=MOST, `1`=unified RSL | `0` |
| `lambdars` | real | RSL influence factor | - |
| `pai_opt` | integer | PAI method: `0`=fixed, `1`=calculated, `2`=from LAI, `3`=user-set | `0` |
| `pai_set` | real | User-set PAI value (used if `pai_opt=3`) | `4.0` |
| `z0_opt` | integer | Surface roughness: `0`=model input, `1`=vegetation-dependent | `0` |

### Wind Adjustment Factor (WAF)

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanwaf` | logical | Enable WAF calculations | `.FALSE.` |
| `dx_opt` | integer | Grid resolution: `0`=calculate from lon, `1`=user-set | `0` |
| `dx_set` | real | Grid resolution (m) | - |
| `flameh_opt` | integer | Flame height method (0-5) | - |
| `flameh_cal` | integer | FRP calculation: `0`=Table 1 method, `1`=Table 2 method | - |
| `flameh_set` | real | Flame height (m) or fraction of canopy height | - |
| `frp_fac` | real | FRP tuning factor | `1.0` |

**Flame Height Options (`flameh_opt`):**
- `0`: Calculate from FRP (fire intensity)
- `1`: User-set flame height
- `2`: FRP where available, user-set elsewhere
- `3`: Fraction of canopy height override
- `4`: FRP where available, fraction override elsewhere
- `5`: FRP/intensity dependent with crown fire detection

### Vertical Diffusion

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcaneddy` | logical | Enable eddy diffusivity calculations | `.FALSE.` |

### Photolysis

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanphot` | logical | Enable photolysis attenuation | `.FALSE.` |

## Biogenic Emissions Options

### Basic Settings

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanbio` | logical | Enable biogenic emissions | `.FALSE.` |
| `bio_cce` | real | Canopy environment coefficient | `0.21` |
| `biospec_opt` | integer | Species output: `0`=all, `1-19`=single species | `0` |
| `biovert_opt` | integer | Vertical summing: `0`=none, `1`=LAD-weighted, `2`=Gaussian, `3`=even | `0` |

### Loss and Chemistry Factors

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `loss_opt` | integer | Canopy loss factor: `0`=off, `1`=calculated, `2`=constant | `0` |
| `loss_set` | real | Constant loss factor (used if `loss_opt=2`) | `0.96` |
| `loss_ind` | integer | Species for loss: `0`=all, `>0`=specific species ID | `0` |
| `lifetime` | real | Chemical lifetime (s) for loss calculations | `3600` |
| `co2_opt` | integer | CO₂ inhibition: `0`=Possell & Hewitt, `1`=Wilkinson, `2`=off | `0` |
| `co2_set` | real | CO₂ concentration (ppmv) | - |

### Environmental Response Factors

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `leafage_opt` | integer | Leaf age response: `0`=on, `1`=off | `1` |
| `lai_tstep` | integer | LAI data interval (s) - e.g., 86400 for daily | - |
| `hist_opt` | integer | Historical averaging: `0`=off, `1`=on | `0` |
| `soim_opt` | integer | Soil moisture response: `0`=on, `1`=off | `1` |
| `soild1` | real | Soil layer 1 depth (cm) | `5.0` |
| `soild2` | real | Soil layer 2 depth (cm) | `25.0` |
| `soild3` | real | Soil layer 3 depth (cm) | `70.0` |
| `soild4` | real | Soil layer 4 depth (cm) | `150.0` |

### Stress Factors

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `aq_opt` | integer | Air quality stress: `0`=calculated, `1`=user-set, `2`=off | `2` |
| `w126_set` | real | Ozone W126 value (ppm-hours) | - |
| `ht_opt` | integer | High temperature stress: `0`=on, `1`=off | `1` |
| `lt_opt` | integer | Low temperature stress: `0`=on, `1`=off | `1` |
| `hw_opt` | integer | High wind stress: `0`=on, `1`=off | `1` |

## Gas Dry Deposition Options

### Basic Settings

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `ifcanddepgas` | logical | Enable gas dry deposition | `.FALSE.` |
| `ddepspecgas_opt` | integer | Species output: `0`=all, `1-31`=single species | `0` |
| `chemmechgas_opt` | integer | Chemical mechanism: `0`=RACM2 | `0` |
| `chemmechgas_tot` | integer | Total species in mechanism | `31` |

### Physical Settings

| Option | Type | Description | Default |
|--------|------|-------------|---------|
| `hyblev1` | real | Height of first hybrid level (m) | `20.0` |
| `snowc_set` | real | Snow cover threshold (%) | `50.0` |
| `icec_set` | real | Ice cover threshold (%) | `50.0` |
| `gamma_set` | real | Building surface reaction probability | `5.0e-5` |
| `Ramin_set` | real | Minimum aerodynamic resistance (s/m) | `10.0` |

## Output Species Tables

### Biogenic Emissions Species (Table 1)

| ID | Variable | Description | Units |
|----|----------|-------------|--------|
| 1 | `emi_isop` | Isoprene | kg m⁻³ s⁻¹ |
| 2 | `emi_myrc` | Myrcene | kg m⁻³ s⁻¹ |
| 3 | `emi_sabi` | Sabinene | kg m⁻³ s⁻¹ |
| 4 | `emi_limo` | Limonene | kg m⁻³ s⁻¹ |
| 5 | `emi_care` | 3-Carene | kg m⁻³ s⁻¹ |
| 6 | `emi_ocim` | t-beta-Ocimene | kg m⁻³ s⁻¹ |
| 7 | `emi_bpin` | beta-Pinene | kg m⁻³ s⁻¹ |
| 8 | `emi_apin` | alpha-Pinene | kg m⁻³ s⁻¹ |
| 9 | `emi_mono` | Other Monoterpenes (34 compounds) | kg m⁻³ s⁻¹ |
| 10 | `emi_farn` | alpha-Farnesene | kg m⁻³ s⁻¹ |
| 11 | `emi_cary` | beta-Caryophyllene | kg m⁻³ s⁻¹ |
| 12 | `emi_sesq` | Other Sesquiterpenes (30 compounds) | kg m⁻³ s⁻¹ |
| 13 | `emi_mbol` | 2-Methyl-3-buten-2-ol | kg m⁻³ s⁻¹ |
| 14 | `emi_meth` | Methanol | kg m⁻³ s⁻¹ |
| 15 | `emi_acet` | Acetone | kg m⁻³ s⁻¹ |
| 16 | `emi_co` | Carbon Monoxide | kg m⁻³ s⁻¹ |
| 17 | `emi_bvoc` | Bi-Directional VOCs (5 compounds) | kg m⁻³ s⁻¹ |
| 18 | `emi_svoc` | Stress VOCs (15 compounds) | kg m⁻³ s⁻¹ |
| 19 | `emi_ovoc` | Other VOCs (49 compounds) | kg m⁻³ s⁻¹ |

### Gas Dry Deposition Species - RACM2 (Table 2)

| ID | Variable | Description | Units |
|----|----------|-------------|--------|
| 1 | `ddep_so2` | Sulfur Dioxide | cm s⁻¹ |
| 2 | `ddep_no` | Nitric Oxide | cm s⁻¹ |
| 3 | `ddep_no2` | Nitrogen Dioxide | cm s⁻¹ |
| 4 | `ddep_nh3` | Ammonia | cm s⁻¹ |
| 5 | `ddep_n2o5` | Dinitrogen Pentoxide | cm s⁻¹ |
| 6 | `ddep_hno3` | Nitric Acid | cm s⁻¹ |
| 7 | `ddep_hono` | Nitrous Acid | cm s⁻¹ |
| 8 | `ddep_pna` | Peroxynitric Acid | cm s⁻¹ |
| 9 | `ddep_h2o2` | Hydrogen Peroxide | cm s⁻¹ |
| 10 | `ddep_o3` | Ozone | cm s⁻¹ |
| 11 | `ddep_form` | Formaldehyde | cm s⁻¹ |
| 12 | `ddep_ald` | Acetaldehyde | cm s⁻¹ |
| 13 | `ddep_ora2` | Organic Acid | cm s⁻¹ |
| 14 | `ddep_aacd` | Acetic Acid | cm s⁻¹ |
| 15 | `ddep_pacd` | Peroxyacetic Acid | cm s⁻¹ |
| 16 | `ddep_op1` | Methyl Hydroperoxide | cm s⁻¹ |
| 17 | `ddep_op2` | Higher Organic Peroxide | cm s⁻¹ |
| 18 | `ddep_pan` | Peroxyacetyl Nitrate | cm s⁻¹ |
| 19 | `ddep_mpan` | Peroxymethacryloyl Nitrate | cm s⁻¹ |
| 20 | `ddep_ntr` | Organic Nitrate | cm s⁻¹ |
| 21 | `ddep_facd` | Formic Acid | cm s⁻¹ |
| 22 | `ddep_co` | Carbon Monoxide | cm s⁻¹ |
| 23 | `ddep_eth` | Ethane | cm s⁻¹ |
| 24 | `ddep_hc3` | Alkanes | cm s⁻¹ |
| 25 | `ddep_hc5` | Higher Alkanes | cm s⁻¹ |
| 26 | `ddep_hc8` | Higher Alkanes | cm s⁻¹ |
| 27 | `ddep_ete` | Ethene | cm s⁻¹ |
| 28 | `ddep_olt` | Terminal Olefins | cm s⁻¹ |
| 29 | `ddep_oli` | Internal Olefins | cm s⁻¹ |
| 30 | `ddep_tol` | Toluene | cm s⁻¹ |
| 31 | `ddep_xyl` | Xylene | cm s⁻¹ |

## Configuration Tips and Best Practices

### Simulation Duration
- **Short runs (< 24 hours)**: Use `hist_opt = 0`
- **Medium runs (1-10 days)**: Use `hist_opt = 1`
- **Long runs (> 10 days)**: Include 240-hour spin-up period

### Vertical Resolution
- **WAF calculations**: Use `modres ≤ 0.5m`
- **Standard applications**: Use `modres = 0.5-2.0m`
- **Large domains**: Use `modres = 1.0-2.0m`

### Biogenic Emissions
- Always use `bio_cce = 0.21` (recommended)
- Enable `soim_opt = 0` for soil moisture effects
- Use `hist_opt = 1` for runs > 24 hours
- Set appropriate `lai_tstep` for your LAI data frequency

### Dry Deposition
- Use default `gamma_set = 5.0e-5` for urban areas
- Adjust snow/ice thresholds (`snowc_set`, `icec_set`) for your domain
- Set `hyblev1` to match your model's first level height

### Wind Calculations
- Set `href_set = 10.0m` (recommended)
- Use `z0ghc = 0.1` for typical forest conditions
- Avoid `href_set` approaching surface roughness length

## Validation and Testing

### Check Configuration
```bash
# Validate namelist syntax
./canopy --check-namelist

# Test short simulation
# Set ntime = 1 for quick validation
```

### Common Validation Steps
1. Verify all required parameters are set
2. Check that file paths exist and are readable
3. Ensure grid dimensions match input files
4. Validate time settings are consistent
5. Check that enabled components have required parameters

### Parameter Ranges
- `modres`: 0.1 - 5.0 m
- `modlays`: 20 - 500 layers
- `bio_cce`: 0.1 - 0.5
- `href_set`: 5.0 - 50.0 m
- `lai_thresh`: 0.1 - 10.0 m²/m²
- `cf_thresh`: 0.1 - 1.0
- `ch_thresh`: 1.0 - 100.0 m

## Related Documentation

- **[Configuration Guide](../getting-started/configuration.md)** - Basic setup instructions
- **[Examples](examples.md)** - Practical configuration examples
- **[Science Documentation](../science/model-description.md)** - Physics behind the parameters
- **[Troubleshooting](troubleshooting.md)** - Problem-solving guide
