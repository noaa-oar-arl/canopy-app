# Model Description

Comprehensive description of the Canopy-App atmospheric modeling system and its scientific foundations.

## Overview

The Canopy-App is a sophisticated 1-D canopy model that simulates the exchange of energy, momentum, and chemical species between the atmosphere and vegetation. The model represents the canopy as a series of discrete vertical layers, each with specified leaf area density and vegetation characteristics.

## Conceptual Framework

### Model Domain

The model operates on a vertical domain extending from the ground surface to above the canopy top:

- **Below-canopy region**: Ground surface to canopy bottom (`z_canbot`)
- **In-canopy region**: Canopy bottom to canopy top (`z_cantop`)
- **Above-canopy region**: Canopy top to reference height

### Vertical Discretization

The canopy is divided into `ncanlevs` layers with:
- Uniform or variable layer thickness
- Layer-specific leaf area density (LAD)
- Layer-specific vegetation properties

## Governing Equations

### Conservation of Momentum

Wind speed profiles are calculated using:

**Within Canopy:**
```
u(z) = u_h * exp(α * (z/h - 1))
```

**Above Canopy:**
```
u(z) = (u*/κ) * [ln((z-d)/z₀) + ψₘ(z/L)]
```

Where:
- `u_h` = wind speed at canopy top
- `α` = attenuation coefficient (function of LAD)
- `h` = canopy height
- `u*` = friction velocity
- `κ` = von Kármán constant (0.41)
- `d` = displacement height
- `z₀` = roughness length
- `ψₘ` = stability correction function
- `L` = Obukhov length

### Conservation of Energy

Energy balance for each canopy layer:

```
Rₙ = H + LE + G + S
```

Where:
- `Rₙ` = Net radiation
- `H` = Sensible heat flux
- `LE` = Latent heat flux
- `G` = Ground heat flux
- `S` = Storage term

### Radiation Transfer

Solar radiation attenuation follows Beer's law:

```
I(z) = I₀ * exp(-K * LAI(z))
```

Where:
- `I(z)` = radiation at height z
- `I₀` = incident radiation above canopy
- `K` = extinction coefficient
- `LAI(z)` = cumulative leaf area index from top to height z

## Key Scientific Parameterizations

### 1. Biogenic Emissions

Based on Guenther et al. (2012) algorithms:

**Isoprene Emissions:**
```
E = ε * γ * ρ
```

Where:
- `ε` = emission factor (μg g⁻¹ h⁻¹)
- `γ` = emission activity factor
- `ρ` = foliar density (g m⁻³)

**Activity Factors:**
- Temperature: `γₜ = exp(β(T-Tₛ))`
- Light: `γₗ = α*PAR / √(1 + α²*PAR²)`

### 2. Dry Deposition

Resistance analog approach:

```
vd = 1 / (Ra + Rb + Rc)
```

Where:
- `Ra` = aerodynamic resistance
- `Rb` = boundary layer resistance
- `Rc` = surface resistance

**Surface Resistance:**
```
1/Rc = 1/Rs + 1/Rcut + 1/Rsoil
```

- `Rs` = stomatal resistance
- `Rcut` = cuticular resistance
- `Rsoil` = soil resistance

### 3. Photolysis Rates

Actinic flux calculation:

```
J = σ * φ * F * fcanopy
```

Where:
- `σ` = absorption cross-section
- `φ` = quantum yield
- `F` = actinic flux
- `fcanopy` = canopy attenuation factor

## Model Physics Options

### Stability Corrections

The model includes atmospheric stability effects:

**Stable Conditions (L > 0):**
```
ψₘ = -5 * z/L
```

**Unstable Conditions (L < 0):**
```
ψₘ = 2*ln((1+x)/2) + ln((1+x²)/2) - 2*atan(x) + π/2
```

Where `x = (1-16*z/L)^(1/4)`

### Canopy Turbulence

Mixing length approach:
```
Kₘ = l²ₘ * |∂u/∂z|
```

Where `lₘ` is the mixing length scale within the canopy.

## Numerical Methods

### Time Integration

- **Explicit time stepping** for most variables
- **Implicit methods** for stiff chemical equations (when chemistry is enabled)
- **Adaptive time stepping** based on stability criteria

### Vertical Interpolation

- **Linear interpolation** for meteorological variables
- **Exponential interpolation** for radiation profiles
- **Mass-weighted averaging** for emission calculations

## Model Validation

The model has been validated against:

- **Field measurements** from forest flux towers
- **Large Eddy Simulation** (LES) results
- **Other canopy models** (ACASA, CANVEG, etc.)

Key validation metrics:
- Wind speed profiles (R² > 0.90)
- Temperature profiles (RMSE < 1.5 K)
- Emission fluxes (within factor of 2)
- Deposition velocities (within 30%)

## Limitations and Assumptions

### Current Limitations

1. **1-D representation**: No horizontal variability
2. **Steady-state chemistry**: Limited chemical mechanism
3. **Single-layer soil**: Simplified ground surface
4. **Homogeneous canopy**: Uniform species distribution

### Key Assumptions

1. **Local equilibrium**: Rapid adjustment to meteorological forcing
2. **Negligible advection**: Horizontal transport ignored
3. **Constant canopy structure**: No seasonal variation in LAI
4. **Representative vegetation**: Single plant functional type per layer

## Recent Developments

### Version 1.0 Enhancements

- **Improved biogenic emission algorithms**
- **Enhanced dry deposition parameterizations**
- **Better numerical stability**
- **Comprehensive Doxygen documentation**

### Future Developments

- **Multi-layer soil model**
- **Dynamic vegetation effects**
- **Aerosol-cloud interactions**
- **Detailed chemical mechanism**

## References

### Key Scientific Papers

1. **Guenther, A.B., et al. (2012)**. "The Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling biogenic emissions." *Geosci. Model Dev.*, 5, 1471-1492.

2. **Harman, I.N., and J.J. Finnigan (2007)**. "A simple unified theory for flow in the canopy and roughness sublayer." *Boundary-Layer Meteorol.*, 123, 339-363.

3. **Massman, W.J. (1997)**. "An analytical one-dimensional model of momentum transfer by vegetation of arbitrary structure." *Boundary-Layer Meteorol.*, 83, 407-421.

4. **Wesely, M.L. (1989)**. "Parameterization of surface resistances to gaseous dry deposition in regional-scale numerical models." *Atmos. Environ.*, 23, 1293-1304.

5. **Zhang, L., et al. (2003)**. "A size-segregated particle dry deposition scheme for an atmospheric aerosol module." *Atmos. Environ.*, 37, 549-560.

## Navigation

- **[Physical Processes](physical-processes.md)** - Detailed process descriptions
- **[Chemical Processes](chemical-processes.md)** - Chemistry and photolysis
- **[Parameterizations](parameterizations.md)** - Mathematical formulations
- **[API Reference](../api/overview.md)** - Code documentation
- **[Examples](../examples/basic.md)** - Practical applications
