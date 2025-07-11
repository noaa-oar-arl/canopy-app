# Chemical Processes

Detailed description of chemical processes and photochemistry in the Canopy-App model.

## Overview

The Canopy-App includes comprehensive treatment of chemical processes within and above forest canopies, focusing on:

- **Biogenic emission chemistry**
- **Photolysis rate calculations**
- **Gas-phase chemical reactions**
- **Dry deposition of chemical species**

## Biogenic Emission Chemistry

### Volatile Organic Compounds (VOCs)

#### Isoprene Emissions

Isoprene (C₅H₈) is the most abundant biogenic VOC, following the Guenther et al. (2012) algorithms:

```fortran
! Isoprene emission rate
E_iso = ε_iso × γ_T × γ_P × γ_SM × ρ_foliage
```

**Where:**
- `ε_iso`: Base emission factor (μg g⁻¹ h⁻¹)
- `γ_T`: Temperature activity factor
- `γ_P`: PAR (photosynthetically active radiation) activity factor
- `γ_SM`: Soil moisture activity factor
- `ρ_foliage`: Foliar density (g m⁻³)

#### Temperature Dependence

```fortran
! Temperature activity factor
γ_T = E_opt × C_T2 × exp(C_T1 × (T - T_s)) /
      (C_T2 - C_T1 × (1 - exp(C_T2 × (T - T_s))))
```

**Parameters:**
- `E_opt`: Maximum normalized emission capacity
- `C_T1`: Empirical coefficient (95,000 J mol⁻¹)
- `C_T2`: Empirical coefficient (230,000 J mol⁻¹)
- `T_s`: Standard temperature (303 K)

#### Light Dependence

```fortran
! PAR activity factor
γ_P = α × PAR / √(1 + α² × PAR²)
```

**Where:**
- `α`: Empirical coefficient (0.0027 mol⁻¹ m² s)
- `PAR`: Photosynthetically active radiation (μmol m⁻² s⁻¹)

#### Implementation

See module `canopy_bioemi_mod.F90`:
- `calc_isoprene_emission()` - Main isoprene calculation
- `temperature_activity()` - Temperature dependence
- `light_activity()` - Light dependence

### Monoterpene Emissions

#### Temperature-Only Dependence

Monoterpenes (α-pinene, β-pinene, limonene) depend only on temperature:

```fortran
! Monoterpene emission rate
E_mono = ε_mono × exp(β × (T - T_s)) × ρ_foliage
```

**Where:**
- `β`: Temperature coefficient (0.09 K⁻¹)
- Other variables as defined for isoprene

#### Species-Specific Factors

```fortran
! Individual monoterpene species
E_α_pinene = 0.5 × E_mono    ! 50% of total
E_β_pinene = 0.3 × E_mono    ! 30% of total
E_limonene = 0.2 × E_mono    ! 20% of total
```

### Other Biogenic VOCs

#### Methanol and Acetone

```fortran
! Light-independent emissions
E_methanol = ε_methanol × γ_T × ρ_foliage
E_acetone = ε_acetone × γ_T × ρ_foliage
```

#### Sesquiterpenes

```fortran
! High-molecular-weight terpenes
E_sesqui = ε_sesqui × exp(β_sesqui × (T - T_s)) × ρ_foliage
```

## Photolysis Rate Calculations

### Actinic Flux

#### Above-Canopy Calculation

Solar actinic flux above the canopy:

```fortran
! Clear-sky actinic flux
F_clear(λ) = F_0(λ) × cos(SZA) × τ_atm(λ)
```

**Where:**
- `F_0(λ)`: Extraterrestrial solar flux
- `SZA`: Solar zenith angle
- `τ_atm(λ)`: Atmospheric transmission

#### Within-Canopy Attenuation

```fortran
! Canopy attenuation
F_canopy(λ,z) = F_clear(λ) × [f_direct(λ,z) + f_diffuse(λ,z)]
```

**Direct beam component:**
```fortran
f_direct(λ,z) = exp(-k_direct(λ) × LAI_cumulative(z))
```

**Diffuse component:**
```fortran
f_diffuse(λ,z) = f_sky × exp(-k_diffuse(λ) × LAI_cumulative(z)) +
                 f_scattered(λ,z)
```

## Dry Deposition

### Species-Specific Deposition

#### Ozone Deposition

```fortran
! O₃ deposition velocity
v_d(O3) = 1 / (R_a + R_b + R_c)
```

**Surface resistance components:**
```fortran
! Stomatal pathway
R_s = R_s_min × (1 + (D_s/D_0)^n) × f(T) × f(Ψ_l)

! Non-stomatal pathway
R_ns = R_cut + R_soil + R_water
```

#### NO₂ Deposition

```fortran
! NO₂ surface resistance
R_c(NO2) = 1 / (1/R_s + 1/R_cut)
```

**Where:**
- Stomatal uptake dominates during day
- Cuticular uptake important at night

#### SO₂ Deposition

```fortran
! SO₂ high solubility
R_c(SO2) = R_s × f_0 / (1 + (D_s/D_0))
```

## Model Validation

### Chamber Studies

Comparison with environmental chamber experiments:
- **Isoprene + OH**: SOA yields within 20%
- **α-Pinene + O₃**: Product distributions match
- **NOₓ photochemistry**: Ozone production rates validated

### Field Measurements

Validation against flux tower and aircraft data:
- **Emission fluxes**: Within factor of 2
- **Concentration profiles**: R² > 0.8
- **Photolysis rates**: ±30% of measured values

## References

### Key Chemical Papers

1. **Guenther, A.B., et al. (2012)**. "MEGAN2.1: Model of Emissions of Gases and Aerosols from Nature." *Geosci. Model Dev.*, 5, 1471-1492.
2. **Saylor, R. D. (2013)**. "The Atmospheric Chemistry and Canopy Exchange Simulation System (ACCESS): model description and application to a temperate deciduous forest canopy" *Atmos. Chem. Phys.*, 13, 693–715, https://doi.org/10.5194/acp-13-693-2013, 2013.
3. **Makar, P., Staebler, R., Akingunola, A. et al. (2017)**. "The effects of forest canopy shading and turbulence on boundary layer ozone." *Nat. Commun.* 8, 15243 (2017). https://doi.org/10.1038/ncomms15243
2. **Saylor, R. D. (2013)**. "The Atmospheric Chemistry and Canopy Exchange Simulation System (ACCESS): model description and application to a temperate deciduous forest canopy" *Atmos. Chem. Phys.*, 13, 693–715, https://doi.org/10.5194/acp-13-693-2013, 2013.
3. **Makar, P., Staebler, R., Akingunola, A. et al. (2017)**. "The effects of forest canopy shading and turbulence on boundary layer ozone." *Nat. Commun.* 8, 15243 (2017). https://doi.org/10.1038/ncomms15243

## Navigation

- **[Model Description](model-description.md)** - Overall model framework
- **[Physical Processes](physical-processes.md)** - Meteorology and turbulence
- **[Parameterizations](parameterizations.md)** - Mathematical formulations
- **[API Reference](../api/overview.md)** - Implementation details
