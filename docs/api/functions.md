# Function Reference

Comprehensive listing of all functions and subroutines in the Canopy-App system.

!!! info "Auto-Generated Documentation"
    This page provides an index of functions. For detailed documentation with parameters, return values, and examples, see the **[Auto-Generated Doxygen API](../canopy/functions.md)**.

## Quick Access to Function Categories

::: canopy
    selection:
        functions: true
    rendering:
        show_source: false
        show_bases: false
        heading_level: 3

## Function Groups

### Core System Functions
Functions for model initialization, execution, and cleanup.

::: canopy
    selection:
        members:
            - canopy_app
            - MemoryManagement
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4
        show_symbol_type_toc: true

### I/O Functions
File handling and data input/output operations.

::: canopy
    selection:
        members:
            - files_mod
            - ncfio_group
            - txtio_group
            - check_input
    rendering:
        show_source: false
        show_bases: false
        heading_level: 4
        show_symbol_type_toc: true

## Complete Function Listing

For the complete alphabetical listing of all functions with full documentation:

**[📋 Browse All Functions in Doxygen API](../canopy/functions.md)**

## A

### `alloc_canopy_arrays()`
**Module:** `canopy_alloc.F90`
**Purpose:** Allocate memory for all canopy-related arrays
**Parameters:** Integer status return code

### `atmos_stability()`
**Module:** `canopy_canmet_mod.F90`
**Purpose:** Calculate atmospheric stability corrections
**Parameters:** Height, temperature, wind speed profiles

## B

### `bioemi_calculate()`
**Module:** `canopy_bioemi_mod.F90`
**Purpose:** Calculate biogenic emission rates for all species
**Parameters:** Temperature, PAR, LAI profiles

### `bioemi_isoprene()`
**Module:** `canopy_bioemi_mod.F90`
**Purpose:** Calculate isoprene emission rates
**Parameters:** Temperature, PAR, emission factors

### `bioemi_monoterpene()`
**Module:** `canopy_bioemi_mod.F90`
**Purpose:** Calculate monoterpene emission rates
**Parameters:** Temperature, emission factors

## C

### `calc_solar_zenith()`
**Module:** `canopy_rad_mod.F90`
**Purpose:** Calculate solar zenith angle for given date/time/location
**Parameters:** Date, time, latitude, longitude

### `canopy_main()`
**Module:** `canopy_app.F90`
**Purpose:** Main program execution routine
**Parameters:** Command line arguments

### `check_input_files()`
**Module:** `canopy_check_input.F90`
**Purpose:** Validate existence and format of input files
**Parameters:** File paths and format specifications

### `check_physics_params()`
**Module:** `canopy_check_input.F90`
**Purpose:** Validate physics parameter ranges and consistency
**Parameters:** Physics parameter arrays

## D

### `date_to_julian()`
**Module:** `canopy_date_mod.F90`
**Purpose:** Convert calendar date to Julian day number
**Parameters:** Year, month, day integers

### `dealloc_canopy_arrays()`
**Module:** `canopy_dealloc.F90`
**Purpose:** Deallocate all dynamically allocated arrays
**Parameters:** Integer status return code

### `drydep_calculate()`
**Module:** `canopy_drydep_mod.F90`
**Purpose:** Calculate dry deposition velocities and fluxes
**Parameters:** Species properties, meteorological conditions

### `drydep_resistance()`
**Module:** `canopy_drydep_mod.F90`
**Purpose:** Calculate resistance components for dry deposition
**Parameters:** Meteorology, surface properties, species data

## E

### `eddy_diffusivity()`
**Module:** `canopy_eddy_mod.F90`
**Purpose:** Calculate turbulent eddy diffusivity profiles
**Parameters:** Wind profiles, stability parameters

## F

### `file_exists()`
**Module:** `canopy_files_mod.F90`
**Purpose:** Check if specified file exists and is readable
**Parameters:** File path string

### `file_format_detect()`
**Module:** `canopy_files_mod.F90`
**Purpose:** Automatically detect file format (NetCDF vs text)
**Parameters:** File path string

## G

### `get_file_extension()`
**Module:** `canopy_utils_mod.F90`
**Purpose:** Extract file extension from path
**Parameters:** File path string

## I

### `init_canopy_grid()`
**Module:** `canopy_init.F90`
**Purpose:** Initialize vertical canopy grid structure
**Parameters:** Canopy height, number of levels

### `init_physics_constants()`
**Module:** `canopy_init.F90`
**Purpose:** Initialize physical constants and default parameters
**Parameters:** None

### `interpolate_profile()`
**Module:** `canopy_utils_mod.F90`
**Purpose:** Interpolate values between vertical levels
**Parameters:** Input/output heights and values

## J

### `julian_to_date()`
**Module:** `canopy_date_mod.F90`
**Purpose:** Convert Julian day number to calendar date
**Parameters:** Julian day number

## L

### `lai_profile()`
**Module:** `canopy_profile_mod.F90`
**Purpose:** Calculate leaf area index profile through canopy
**Parameters:** Total LAI, canopy height, distribution parameters

## N

### `ncf_read_input()`
**Module:** `canopy_ncf_io_mod.F90`
**Purpose:** Read meteorological data from NetCDF input file
**Parameters:** File path, variable arrays

### `ncf_write_output()`
**Module:** `canopy_ncf_io_mod.F90`
**Purpose:** Write model results to NetCDF output file
**Parameters:** File path, output variable arrays

## P

### `photolysis_calculate()`
**Module:** `canopy_phot_mod.F90`
**Purpose:** Calculate photolysis rates for all species
**Parameters:** Solar flux, cross sections, quantum yields

### `photolysis_actinic()`
**Module:** `canopy_phot_mod.F90`
**Purpose:** Calculate actinic flux attenuation through canopy
**Parameters:** Solar zenith angle, LAI profile

## R

### `radiation_transfer()`
**Module:** `canopy_rad_mod.F90`
**Purpose:** Calculate radiation transfer through canopy layers
**Parameters:** Solar flux, LAI profile, extinction coefficients

### `read_namelist_inputs()`
**Module:** `canopy_readnml.F90`
**Purpose:** Read and parse canopy_inputs namelist
**Parameters:** Namelist file path

### `read_namelist_options()`
**Module:** `canopy_readnml.F90`
**Purpose:** Read and parse canopy_options namelist
**Parameters:** Namelist file path

### `read_point_data()`
**Module:** `canopy_read_txt.F90`
**Purpose:** Read point measurement data from text file
**Parameters:** File path, data arrays

## S

### `stability_function()`
**Module:** `canopy_canmet_mod.F90`
**Purpose:** Calculate stability correction functions
**Parameters:** Stability parameter, height

### `surface_resistance()`
**Module:** `canopy_drydep_mod.F90`
**Purpose:** Calculate surface resistance for gas deposition
**Parameters:** Species properties, environmental conditions

## T

### `temperature_profile()`
**Module:** `canopy_canmet_mod.F90`
**Purpose:** Calculate temperature profile through canopy
**Parameters:** Reference temperature, heat flux, stability

### `txt_read_input()`
**Module:** `canopy_txt_io_mod.F90`
**Purpose:** Read input data from ASCII text files
**Parameters:** File path, format specifications

### `txt_write_output()`
**Module:** `canopy_txt_io_mod.F90`
**Purpose:** Write output data to ASCII text files
**Parameters:** File path, output arrays, format specifications

## U

### `unit_conversion()`
**Module:** `canopy_utils_mod.F90`
**Purpose:** Convert between different units
**Parameters:** Input value, input units, output units

## V

### `validate_date_string()`
**Module:** `canopy_date_mod.F90`
**Purpose:** Validate date string format and values
**Parameters:** Date string

### `vertical_interpolation()`
**Module:** `canopy_coord_mod.F90`
**Purpose:** Interpolate between vertical coordinate systems
**Parameters:** Input/output coordinate arrays

## W

### `wind_profile()`
**Module:** `canopy_canmet_mod.F90`
**Purpose:** Calculate wind speed profile through canopy
**Parameters:** Reference wind speed, canopy structure, stability

### `write_point_output()`
**Module:** `canopy_write_txt.F90`
**Purpose:** Write point output data to text file
**Parameters:** File path, point data arrays

## Usage Notes

### Function Call Patterns
Most Canopy-App subroutines follow this pattern:
```fortran
call subroutine_name(input_params, output_params, status)
if (status /= 0) then
    ! Handle error condition
    write(*,*) 'Error in subroutine_name: ', status
endif
```

### Error Handling
- Return status codes: 0 = success, non-zero = error
- Error messages written to standard output
- Critical errors may halt program execution

### Parameter Conventions
- `intent(in)` - Input parameters
- `intent(out)` - Output parameters
- `intent(inout)` - Input/output parameters
- `optional` - Optional parameters

## Next Steps

- 📋 [Browse modules](modules.md)
- 📊 [View variables](variables.md)
- 🔗 [Complete API documentation](../canopy/links.md)
- 💻 [See usage examples](../examples/basic.md)
