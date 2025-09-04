# Variables Reference

Complete reference for all variables, parameters, and constants in the Canopy-App system.

!!! info "Auto-Generated Documentation"
    This page provides an overview of variables. For detailed documentation with types, units, and usage, see the **[Auto-Generated Doxygen Variables](../canopy/variables.md)**.

## Quick Access to Variable Categories

::: canopy
    selection:
        variables: true
    rendering:
        show_source: false
        show_bases: false
        heading_level: 3

## Variable Groups

### Physical Constants
::: canopy
    selection:
        members:
            - const_mod
            - fundamental_const
            - atmos_const
            - bioemi_const
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4

Physical and mathematical constants used throughout the model.

### Model Configuration Variables
::: canopy
    selection:
        members:
            - canopy_options
            - canopy_variables
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4

Model configuration options and main canopy variables.

### Input/Output Variables
::: canopy
    selection:
        members:
            - InputVariables
            - file_arrays
            - file_paths
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4

Input data arrays and file path variables.

### Canopy Structure Variables
::: canopy
    selection:
        members:
            - CanopyDistribution
            - coord_mod
            - spatial_domain
            - vertical_domain
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4

Variables defining canopy structure and spatial domains.

### Meteorological Variables
::: canopy
    selection:
        members:
            - canopy_met_input
            - CanopyWind
            - CanopyDiffusivity
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4

Meteorological input and calculated variables.

### Process-Specific Variables
::: canopy
    selection:
        members:
            - CanopyBiogenicEmissions
            - CanopyPhotolysis
            - CanopyDryDeposition
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4

Variables for biogenic emissions, photolysis, and dry deposition processes.

## Complete Variable Listing

For the complete listing of all variables with detailed documentation:

**[📋 Browse All Variables in Doxygen API](../canopy/variables.md)**

## Variable Categories by Type

| Category | Key Variables | Description |
|----------|---------------|-------------|
| **Constants** | `pi`, `grav`, `rgas`, `avogad` | Mathematical and physical constants |
| **Configuration** | `opt_canmet`, `opt_bioem`, `opt_drydep` | Model physics options |
| **Dimensions** | `ncanlevs`, `z_cantop`, `z_canbot` | Canopy structure parameters |
| **Meteorology** | `temp_k`, `wspd_ms`, `rh_pct` | Weather variables |
| **Emissions** | `efiso`, `efmono`, `biomass_burn` | Emission factors and rates |
| **Deposition** | `vd_o3`, `vd_no2`, `stomatal_res` | Deposition velocities and resistances |

## Navigation

- **[Module Reference](modules.md)** - Browse by module organization
- **[Function Reference](functions.md)** - Alphabetical function list
- **[Complete Doxygen API](../canopy/links.md)** - Full auto-generated documentation
- **[Usage Examples](../examples/basic.md)** - Variable usage examples
