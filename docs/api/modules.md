# Module Reference

Complete reference for all Fortran modules in the Canopy-App system. This page provides an overview and links to the detailed auto-generated API documentation.

!!! info "Auto-Generated Documentation"
    For detailed API documentation with function signatures, parameters, and examples, see the **[Auto-Generated Doxygen API](../canopy/links.md)**.

## Core System Modules

### Main Application
::: canopy
    selection:
        members:
            - canopy_app
    rendering:
        show_source: true
        show_bases: false
        heading_level: 4

The main program entry point that coordinates all model execution phases.

### Memory Management
::: canopy
    selection:
        members:
            - MemoryManagement
    rendering:
        show_source: true
        show_bases: false
        heading_level: 4

Dynamic memory allocation and deallocation for all model arrays.

Reads and parses Fortran namelist configuration files.

**Supported Namelists:**
- `&canopy_inputs` - Input/output file specifications
- `&canopy_options` - Physics options and switches
- `&canopy_physics` - Physical parameter settings

---

### canopy_files_mod.F90
**File Management**

Utilities for file path handling, existence checking, and format detection.

**Key Functions:**
- File existence validation
- Path resolution and cleanup
- Format detection (NetCDF vs text)

---

### canopy_ncf_io_mod.F90
**NetCDF Input/Output**

Complete NetCDF file handling for both input and output operations.

**Capabilities:**
- Read meteorological input data
- Write model output in CF-compliant format
- Handle time series and profile data
- Error checking and validation

---

### canopy_txt_io_mod.F90
**Text File Input/Output**

ASCII text file I/O operations for simple data exchange.

**Features:**
- Column-based data reading
- Formatted output generation
- Header processing
- Quality control checks

## Physical Process Modules

### canopy_canmet_mod.F90
**Canopy Meteorology**

Computes meteorological variables within and above the canopy.

**Calculations:**
- Wind speed profiles using similarity theory
- Temperature and humidity gradients
- Atmospheric stability corrections
- Turbulent mixing parameterizations

**Key Variables:**
- `wind_prof(nlev)` - Wind speed profile
- `temp_prof(nlev)` - Temperature profile
- `rh_prof(nlev)` - Relative humidity profile

---

### canopy_rad_mod.F90
**Radiation Transfer**

Solar and longwave radiation calculations through the canopy.

**Capabilities:**
- Beer's law attenuation
- Direct and diffuse radiation separation
- Solar zenith angle calculations
- Photosynthetically active radiation (PAR)

**Key Variables:**
- `solar_flux(nlev)` - Solar radiation by level
- `par_flux(nlev)` - PAR availability
- `zenith_angle` - Solar zenith angle

---

### canopy_bioemi_mod.F90
**Biogenic Emissions**

Vegetation emissions of volatile organic compounds.

**Emission Types:**
- Temperature-dependent emissions (monoterpenes)
- Light-dependent emissions (isoprene)
- Stress-induced emissions
- Species-specific emission factors

**Key Variables:**
- `emis_isop(nlev)` - Isoprene emission rates
- `emis_mono(nlev)` - Monoterpene emission rates
- `emis_total` - Total canopy emissions

---

### canopy_drydep_mod.F90
**Dry Deposition**

Removal of atmospheric gases and particles by vegetation.

**Methods:**
- Multi-layer resistance networks
- Big-leaf approximations
- Stomatal and cuticular pathways
- Environmental controls

**Key Variables:**
- `vdep_gas(nspec)` - Deposition velocities by species
- `flux_dep(nspec)` - Deposition fluxes
- `resist_stom(nlev)` - Stomatal resistance

---

### canopy_phot_mod.F90
**Photolysis Rates**

Photolysis rate calculations for atmospheric chemistry.

**Features:**
- Actinic flux attenuation
- Species-specific cross sections
- Canopy shading effects
- Wavelength-dependent calculations

**Key Variables:**
- `jvals(nphot,nlev)` - Photolysis rates by species and level
- `actinic_flux(nlev)` - Actinic flux profile

## Supporting Modules

### canopy_const_mod.F90
**Physical Constants**

Fundamental physical and mathematical constants used throughout the model.

**Categories:**
- Physical constants (gravity, gas constant, etc.)
- Conversion factors
- Mathematical constants
- Default parameter values

---

### canopy_utils_mod.F90
**Utility Functions**

General-purpose utility functions and mathematical operations.

**Functions:**
- Array operations and statistics
- Unit conversions
- Mathematical utilities
- String processing

---

### canopy_coord_mod.F90
**Coordinate Systems**

Coordinate system transformations and grid utilities.

**Features:**
- Vertical coordinate transformations
- Grid spacing calculations
- Level interpolation routines

---

### canopy_date_mod.F90
**Date and Time Handling**

Date/time parsing, validation, and calculations.

**Capabilities:**
- Date string parsing
- Julian day calculations
- Time zone handling
- Solar angle computations

## Variable and Parameter Modules

### canopy_canvars_mod.F90
**Canopy Variables**

Central repository for all canopy-related variables and arrays.

### canopy_canopts_mod.F90
**Canopy Options**

Model configuration options and physics switches.

### canopy_bioparm_mod.F90
**Biogenic Parameters**

Parameters for biogenic emission calculations.

## Next Steps

- 🔍 [View function index](functions.md)
- 📊 [Browse variable reference](variables.md)
- 🔗 [Access complete Doxygen API](../canopy/links.md)
- 📖 [Read user guide](../user-guide/overview.md)
