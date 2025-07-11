# Physical Processes

Detailed description of the physical processes simulated in the Canopy-App model.

## Canopy Meteorology

### Wind Flow Through Vegetation

#### Above-Canopy Flow

The wind profile above the canopy follows Monin-Obukhov similarity theory:

```fortran
! Above-canopy wind profile
u(z) = (u_star / kappa) * [log((z - d) / z0) + psi_m(z/L)]
```

**Parameters:**
- `u_star`: Friction velocity (m/s)
- `kappa`: von Kármán constant (0.41)
- `d`: Displacement height (m)
- `z0`: Roughness length (m)
- `L`: Obukhov length (m)
- `psi_m`: Momentum stability function

#### Within-Canopy Flow

Wind speed decreases exponentially through the canopy:

```fortran
! In-canopy wind profile
u(z) = u_h * exp(alpha * (z/h - 1.0))
```

**Where:**
- `alpha = LAD_profile * Cd * h` (attenuation coefficient)
- `LAD_profile`: Leaf area density profile (m²/m³)
- `Cd`: Drag coefficient (~0.2-0.3)
- `h`: Canopy height (m)

#### Implementation

See module `canopy_wind_mod.F90`:
- `calc_wind_most()` - Main wind calculation routine

### Temperature Profiles

#### Energy Balance

<!--#### Energy Balance

```fortran
! Layer energy balance
Rn(z) = H(z) + LE(z) + storage_term(z)
```-->
- `Rn`: Net radiation (W/m²)
<!--**Components:**
- `LE`: Latent heat flux (W/m²) -->

#### Sensible Heat Flux

- `LE`: Latent heat flux (W/m²) -->

<!--#### Sensible Heat Flux
! Sensible heat flux
H(z) = -rho * cp * Kh(z) * dT/dz
```

**Where:**
- `Kh`: Eddy diffusivity for heat (m²/s)
- `dT/dz`: Temperature gradient (K/m) -->
- `rho`: Air density (kg/m³)
<!-- - `cp`: Specific heat of air (J/kg/K) -->

<!-- #### Implementation

See module `canopy_canmet_mod.F90`:
- `calc_temperature_profile()` - Temperature calculations
- `calc_heat_flux()` - Sensible heat flux
- `calc_eddy_diffusivity()` - Turbulent mixing -->

<!-- ### Humidity and Latent Heat

#### Evapotranspiration

Total latent heat flux includes:

1. **Transpiration** from stomata
2. **Evaporation** from wet surfaces
3. **Soil evaporation** from ground

```fortran
! Total latent heat flux
LE_total = LE_transpiration + LE_evaporation + LE_soil
``` -->

<!-- #### Stomatal Conductance

Based on environmental controls:

```fortran
! Stomatal conductance
gs = gs_max * f_light * f_temp * f_humidity * f_co2
```

**Environmental factors:**
- `f_light`: Light response function
- `f_temp`: Temperature response function
- `f_humidity`: Humidity stress function
- `f_co2`: CO₂ response function -->

## Radiation Transfer

<!-- ### Solar Radiation Components -->

<!-- #### Direct and Diffuse Radiation

Solar radiation is separated into:
- **Direct beam radiation**: `I_direct`
- **Diffuse radiation**: `I_diffuse`
- **Scattered radiation**: `I_scattered` -->

#### Photosynthetically Active Radiation (PAR)

PAR calculation for photosynthesis and biogenic emissions:

```fortran
! PAR attenuation through canopy
PAR(z) = PAR_top * exp(-K_par * LAI_cumulative(z))
```

**Where:**
- `K_par`: PAR extinction coefficient
- `LAI_cumulative`: Cumulative LAI from canopy top

#### Implementation

See module `canopy_rad_mod.F90`:
<!-- - `canopy_fsun_clu()` - Main radiation routine -->
- `calc_par_profile()` - PAR calculations
- `calc_extinction_coeff()` - Light extinction
- `canopy_fsun_clu()` - Main radiation routine
<!-- - `canopy_ppfd_exp()` - PPFD calculations -->

<!--```fortran
! Net longwave radiation
Rn_lw = Rn_lw_down - Rn_lw_up
``` -->

<!-- **Components:**
- Atmospheric longwave down
- Canopy longwave emission up
- Multiple scattering within canopy -->

<!-- =#### Sky View Factor

Calculated for each canopy layer:

```fortran
! Sky view factor
svf(z) = exp(-K_lw * LAI_above(z))
``` -->

## Turbulent Transport

### Mixing Length Theory

Eddy diffusivity calculation:

```fortran
! Mixing length approach
Km(z) = l_mix(z)^2 * |du/dz|
```

**Mixing length scale:**
- Above canopy: `l_mix = kappa * z`
- Within canopy: `l_mix = alpha * LAD(z)^(-1)`

#### Stability Corrections

**Stable conditions (L > 0):**
```fortran
psi_m = -5.0 * z / L
psi_h = -5.0 * z / L
```

**Unstable conditions (L < 0):**
```fortran
x = (1.0 - 16.0 * z / L)**0.25
psi_m = 2.0*log((1.0+x)/2.0) + log((1.0+x*x)/2.0) - 2.0*atan(x) + pi/2.0
psi_h = 2.0*log((1.0+x*x)/2.0)
```

#### Implementation

<!-- See module `canopy_eddy_mod.F90`:
- `calc_eddy_diffusivity()` - Main turbulence routine
- `calc_mixing_length()` - Mixing length calculation
- `stability_functions()` - Stability corrections -->

## Boundary Layer Processes
- `canopy_eddyx()` - Main turbulence routine

#### Roughness Parameters

Calculated from canopy structure:

```fortran
! Roughness length
z0 = 0.1 * h

! Displacement height
d = 0.7 * h
```

**Where h is canopy height.**

<!-- #### Heat and Moisture Fluxes

Surface layer fluxes using bulk transfer:

```fortran
! Sensible heat flux
H = rho * cp * CH * U * (Ts - Ta)

! Latent heat flux
LE = rho * lv * CE * U * (qs - qa)
```

**Transfer coefficients:**
- `CH`: Heat transfer coefficient
- `CE`: Moisture transfer coefficient
- `U`: Wind speed
- `Ts, Ta`: Surface and air temperature
- `qs, qa`: Surface and air specific humidity -->

### Canopy-Atmosphere Coupling

#### Feedback Mechanisms

1. **Canopy modification** of atmospheric profiles
2. **Surface flux** responses to atmospheric forcing
3. **Momentum absorption** by vegetation
4. **Scalar transport** through canopy layers

#### Implementation

See modules:
- `canopy_profile_mod.F90` - Vertical profile calculations
- `canopy_utils_mod.F90` - Interpolation and utilities

<!-- ## Soil-Canopy Interactions

### Ground Surface Processes

#### Soil Heat Flux

Simple ground heat flux model:

```fortran
! Soil heat flux
G = -k_soil * dT_soil/dz
```

**Where:**
- `k_soil`: Soil thermal conductivity
- `dT_soil/dz`: Soil temperature gradient

#### Soil Evaporation

Resistance-based approach:

```fortran
! Soil evaporation
E_soil = (e_sat - e_air) / (r_soil + r_aero)
```

**Resistances:**
- `r_soil`: Soil surface resistance
- `r_aero`: Aerodynamic resistance to soil -->

<!-- ### Root Zone Processes

#### Water Uptake

Simplified root water uptake:

```fortran
! Root water uptake
S_root = alpha_root * LAI * (theta - theta_wp)
```

**Where:**
- `alpha_root`: Root efficiency parameter
- `theta`: Soil moisture content
- `theta_wp`: Wilting point -->

## Model Numerics

### Vertical Grid

#### Layer Structure

- **Exponential spacing** for fine resolution near ground
- **Uniform spacing** through main canopy
- **Stretched grid** above canopy

#### Grid Generation

```fortran
! Exponential grid near surface
do k = 1, nlevs_surface
   z(k) = z_min * exp((k-1) * dz_factor)
end do

! Uniform grid in canopy
do k = nlevs_surface+1, nlevs_canopy
   z(k) = z_canbot + (k - nlevs_surface) * dz_uniform
end do
```

### Time Integration

<!-- #### Explicit Schemes

Most variables use explicit time stepping:

```fortran
! Forward Euler
var_new = var_old + dt * tendency
``` -->

#### Stability Criteria

Time step limited by:
- **CFL condition**: `dt < dz / u_max`
- **Diffusion condition**: `dt < dz² / (2 * K_max)`

## Validation Data

### Field Observations

Model validated against measurements from:
- **AmeriFlux tower sites**
- **FLUXNET database**
- **Specialized canopy experiments**

### Typical Performance

| Variable | Correlation (R²) | RMSE |
|----------|-----------------|------|
| Wind Speed | 0.85-0.95 | 0.5-1.0 m/s |
| Temperature | 0.90-0.98 | 1.0-2.0 K |
| Humidity | 0.75-0.90 | 5-10% |
| Radiation | 0.95-0.99 | 10-20 W/m² |

## Navigation

- **[Model Description](model-description.md)** - Overall model framework
- **[Chemical Processes](chemical-processes.md)** - Chemistry and emissions
- **[Parameterizations](parameterizations.md)** - Mathematical details
- **[API Reference](../api/overview.md)** - Implementation details
