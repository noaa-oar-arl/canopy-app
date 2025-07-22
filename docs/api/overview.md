# API Reference Overview

This section provides comprehensive API documentation for all Canopy-App modules, functions, and data structures.

## Documentation Structure

The API documentation is automatically generated from the Fortran source code using Doxygen and organized into the following sections:

### 📋 Quick Navigation
- **[Module List](modules.md)** - All Fortran modules in the codebase
- **[Function Index](functions.md)** - Alphabetical list of all functions and subroutines
- **[Variable Reference](variables.md)** - Global variables and parameters
- **[Generated Doxygen API](../canopy/links.md)** - Complete auto-generated documentation

## Module Categories

### Core System Modules
- **`canopy_app`** - Main program entry point
- **`canopy_init`** - Model initialization routines
- **`canopy_alloc`** - Memory allocation management
- **`canopy_dealloc`** - Memory cleanup routines

### Configuration and I/O
- **`canopy_readnml`** - Namelist reading and parsing
- **`canopy_files_mod`** - File management utilities
- **`canopy_ncf_io_mod`** - NetCDF input/output operations
- **`canopy_txt_io_mod`** - Text file I/O operations
- **`canopy_check_input`** - Input validation and quality control

### Physical Process Modules
- **`canopy_canmet_mod`** - Canopy meteorology calculations
- **`canopy_rad_mod`** - Radiation transfer computations
- **`canopy_bioemi_mod`** - Biogenic emission algorithms
- **`canopy_drydep_mod`** - Dry deposition calculations
- **`canopy_phot_mod`** - Photolysis rate computations

### Supporting Modules
- **`canopy_const_mod`** - Physical and mathematical constants
- **`canopy_utils_mod`** - General utility functions
- **`canopy_coord_mod`** - Coordinate system utilities
- **`canopy_date_mod`** - Date and time handling

## Using the API Documentation

### Finding Functions
1. Browse by **module** to see related functions grouped together
2. Use the **function index** for alphabetical lookup
3. Search the **Doxygen API** for detailed parameter information

### Understanding Parameters
Each function/subroutine is documented with:
- **Purpose** - What the routine does
- **Input parameters** - Required inputs with types and units
- **Output parameters** - Results and their formats
- **Side effects** - Any global state changes
- **Dependencies** - Required modules or external libraries

### Code Examples
Many API entries include usage examples showing:
```fortran
! Example function call
call canopy_function(input_param, output_result, status)
if (status /= 0) then
    write(*,*) 'Error in canopy_function'
endif
```

## API Conventions

### Naming Conventions
- **Modules**: `canopy_<purpose>_mod`
- **Subroutines**: `<action>_<object>` (e.g., `calc_windspeed`)
- **Functions**: `<object>_<property>` (e.g., `temperature_profile`)
- **Variables**: `<type>_<name>` (e.g., `real_windspeed`)

### Parameter Conventions
- **`intent(in)`** - Input parameters (read-only)
- **`intent(out)`** - Output parameters (write-only)
- **`intent(inout)`** - Input/output parameters (read-write)
- **`optional`** - Optional parameters with defaults

### Error Handling
Most routines follow consistent error handling:
- Return status codes (0 = success, non-zero = error)
- Error messages written to standard output
- Graceful degradation when possible

## Common Usage Patterns

### Model Initialization
```fortran
! Standard initialization sequence
call canopy_alloc_arrays(status)
call canopy_read_namelist('namelist.canopy', status)
call canopy_init_physics(status)
call canopy_check_inputs(status)
```

### Main Calculation Loop
```fortran
! Typical calculation sequence
call canopy_calc_meteorology(status)
call canopy_calc_radiation(status)
call canopy_calc_biogenic_emissions(status)
call canopy_calc_dry_deposition(status)
call canopy_write_outputs(status)
```

### Memory Management
```fortran
! Proper cleanup sequence
call canopy_write_final_outputs(status)
call canopy_dealloc_arrays(status)
```

## Performance Considerations

### Computational Efficiency
- Most calculations are vectorized for array operations
- Optional physics modules can be disabled to improve performance
- Memory allocation is optimized for typical canopy layer counts

### Memory Usage
- Arrays are allocated dynamically based on canopy configuration
- Large arrays use double precision only when necessary
- Intermediate arrays are deallocated when no longer needed

## Extending the API

### Adding New Modules
1. Follow the naming convention: `canopy_<name>_mod.F90`
2. Include proper Doxygen documentation headers
3. Implement consistent error handling
4. Add appropriate unit tests

### Modifying Existing Functions
1. Maintain backward compatibility when possible
2. Update documentation for any parameter changes
3. Add version notes for significant modifications
4. Test thoroughly before committing changes

## Next Steps

- 📖 [Browse all modules](modules.md)
- 🔍 [Search functions](functions.md)
- 📊 [View variable reference](variables.md)
- 🔗 [Access complete Doxygen API](../canopy/links.md)
- 💻 [See practical examples](../examples/basic.md)
