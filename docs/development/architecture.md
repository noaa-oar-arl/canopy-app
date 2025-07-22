# Model Architecture

This document describes the software architecture and design principles of the Canopy-App model, providing insights for developers working with or extending the codebase.

## Overall Architecture

### Modular Design

The Canopy-App follows a modular architecture with clear separation of concerns:

```
┌─────────────────┐
│   Main Driver   │  ← canopy_app.F90
└─────────────────┘
         │
    ┌────┴────┐
    │ Control │  ← Initialization, timestepping, finalization
    └────┬────┘
         │
┌────────┴────────┐
│  Core Modules   │  ← Physics, chemistry, I/O modules
└─────────────────┘
         │
┌────────┴────────┐
│ Utility Modules │  ← Constants, utilities, coordinate systems
└─────────────────┘
```

### Module Hierarchy

#### Level 1: Utility Modules
- `canopy_const_mod`: Physical and mathematical constants
- `canopy_coord_mod`: Coordinate system definitions
- `canopy_utils_mod`: Utility functions and procedures
- `canopy_date_mod`: Date/time handling

#### Level 2: Data Modules
- `canopy_canvars_mod`: Canopy variable definitions
- `canopy_canopts_mod`: Model options and configurations
- `canopy_bioparm_mod`: Biophysical parameters

#### Level 3: I/O Modules
- `canopy_files_mod`: File handling
- `canopy_ncf_io_mod`: NetCDF I/O operations
- `canopy_txt_io_mod`: Text file I/O operations

#### Level 4: Physics Modules
- `canopy_profile_mod`: Vertical profile calculations
- `canopy_wind_mod`: Wind profile and turbulence
- `canopy_rad_mod`: Radiation transfer
- `canopy_tleaf_mod`: Leaf temperature calculations
- `canopy_phot_mod`: Photosynthesis calculations

#### Level 5: Chemistry Modules
- `canopy_bioemi_mod`: Biogenic emissions
- `canopy_drydep_mod`: Dry deposition
- `canopy_canmet_mod`: Canopy meteorology

#### Level 6: Main Application
- `canopy_app.F90`: Main program driver

## Design Principles

### 1. Modularity
Each module has a specific responsibility and well-defined interfaces:

```fortran
module canopy_example_mod
  use canopy_const_mod, only: r8  ! Explicit imports
  implicit none
  private  ! Default private scope

  ! Public interface
  public :: example_init, example_calc, example_finalize

  ! Private module variables
  real(r8), allocatable :: private_array(:)

contains

  subroutine example_init(nx, ny)
    integer, intent(in) :: nx, ny
    allocate(private_array(nx*ny))
  end subroutine

end module
```

### 2. Data Encapsulation
Model state variables are encapsulated in derived types:

```fortran
type :: canopy_state_type
  ! Meteorological variables
  real(r8), allocatable :: temp(:,:,:)    ! Temperature
  real(r8), allocatable :: qv(:,:,:)      ! Water vapor
  real(r8), allocatable :: wind_u(:,:,:)  ! U-wind component
  real(r8), allocatable :: wind_v(:,:,:)  ! V-wind component

  ! Canopy variables
  real(r8), allocatable :: lai(:,:)       ! Leaf area index
  real(r8), allocatable :: canht(:,:)     ! Canopy height
end type
```

### 3. Error Handling
Consistent error handling throughout the codebase:

```fortran
subroutine safe_allocate_example(array, n, status, errmsg)
  real(r8), allocatable, intent(out) :: array(:)
  integer, intent(in) :: n
  integer, intent(out) :: status
  character(len=*), intent(out) :: errmsg

  allocate(array(n), stat=status)
  if (status /= 0) then
    errmsg = 'Failed to allocate array in safe_allocate_example'
    return
  end if

  errmsg = ''
end subroutine
```

### 4. Performance Optimization
Code is structured for computational efficiency:

- **Memory locality**: Arrays organized for cache-friendly access patterns
- **Vectorization**: Inner loops designed for compiler vectorization
- **Minimal allocations**: Reuse of temporary arrays where possible

## Memory Management

### Allocation Strategy

The model uses a hierarchical allocation strategy:

1. **Global arrays**: Allocated once during initialization
2. **Temporary arrays**: Allocated/deallocated as needed
3. **Local arrays**: Stack allocation for small, temporary data

### Memory Layout

```fortran
! 3D arrays organized as (x, y, z) for cache efficiency
real(r8), allocatable :: temp(:,:,:)    ! (nlon, nlat, nlev)
real(r8), allocatable :: qv(:,:,:)      ! (nlon, nlat, nlev)

! 2D surface arrays
real(r8), allocatable :: lai(:,:)       ! (nlon, nlat)
real(r8), allocatable :: canht(:,:)     ! (nlon, nlat)
```

### Deallocation
Systematic cleanup in finalization routines:

```fortran
subroutine canopy_dealloc()
  use canopy_canvars_mod
  implicit none

  if (allocated(temp)) deallocate(temp)
  if (allocated(qv)) deallocate(qv)
  ! ... other deallocations

end subroutine
```

## Computational Flow

### Main Program Structure

```fortran
program canopy_app
  ! 1. Initialization phase
  call canopy_init()

  ! 2. Time integration loop
  do timestep = 1, ntimesteps
    call canopy_calcs(timestep)
  end do

  ! 3. Finalization phase
  call canopy_dealloc()
end program
```

### Time Integration

The model uses operator splitting for different processes:

```
Time Step n → n+1:
├── 1. Update meteorology
├── 2. Calculate radiation
├── 3. Compute photosynthesis
├── 4. Update energy balance
├── 5. Calculate turbulence
├── 6. Process chemistry (if enabled)
└── 7. Output results
```

### Operator Splitting

Different physical processes are solved sequentially:

```fortran
subroutine canopy_calcs(itime)
  integer, intent(in) :: itime

  ! Meteorological processes
  call canopy_profile_update()

  ! Radiation calculations
  call canopy_rad_calc()

  ! Leaf-level processes
  call canopy_tleaf_calc()
  call canopy_phot_calc()

  ! Turbulence calculations
  call canopy_eddy_calc()

  ! Chemistry (optional)
  if (chemistry_enabled) then
    call canopy_bioemi_calc()
    call canopy_drydep_calc()
  end if

  ! Output
  call canopy_output(itime)

end subroutine
```

## I/O Architecture

### File Handling Strategy

```
Input Files:
├── Configuration: namelist.canopy
├── Meteorology: *.nc or *.txt
└── Parameters: built-in or external

Processing:
├── Validation
├── Unit conversion
└── Interpolation (if needed)

Output Files:
├── Primary: canopy_output.nc
├── Diagnostics: canopy_diag.nc
└── Text: point_file_*.txt
```

### NetCDF Implementation

The model uses a standardized NetCDF interface:

```fortran
module canopy_ncf_io_mod
  use netcdf
  implicit none

  type :: ncf_file_type
    integer :: ncid
    character(len=256) :: filename
    logical :: is_open
  end type

contains

  subroutine ncf_open_file(file, filename, mode)
    type(ncf_file_type), intent(out) :: file
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: mode

    ! Implementation
  end subroutine

end module
```

## Parallelization Strategy

### Current Implementation
- **OpenMP**: Shared-memory parallelization
- **Thread-safe**: Critical sections protected
- **Load balancing**: Work distributed across threads

### Threading Model

```fortran
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(temp,qv,wind_u,wind_v)
do k = 1, nlev
  do j = 1, nlat
    do i = 1, nlon
      ! Computation for grid point (i,j,k)
      call point_calculation(i, j, k, temp, qv, wind_u, wind_v)
    end do
  end do
end do
!$OMP END PARALLEL DO
```

### Future Scalability
Architecture designed for future enhancements:

- **MPI**: Domain decomposition for distributed memory
- **GPU**: Compute kernels suitable for GPU acceleration
- **Hybrid**: MPI+OpenMP for large-scale systems

## Testing Architecture

### Unit Testing Framework

```fortran
module test_canopy_module
  use canopy_test_framework
  implicit none

contains

  subroutine test_photosynthesis()
    real(r8) :: ppfd, temp, co2, expected, actual

    ! Setup test conditions
    ppfd = 1000.0_r8
    temp = 298.15_r8
    co2 = 400.0_r8
    expected = 25.5_r8

    ! Call function under test
    actual = calc_photosynthesis(ppfd, temp, co2)

    ! Assert result
    call assert_real_equal(expected, actual, 0.1_r8, 'Photosynthesis test')
  end subroutine

end module
```

### Integration Testing

Full model runs with known inputs/outputs for regression testing.

## Configuration Management

### Namelist Structure

```fortran
namelist /canopy_options/ &
  nlat, nlon, nlev, ntime, &
  dx, dy, dz_top, &
  infmt_opt, outfmt_opt, &
  chemistry_enabled, &
  output_freq

namelist /physics_options/ &
  turbulence_scheme, &
  radiation_scheme, &
  photosynthesis_model

namelist /chemistry_options/ &
  biogenic_emissions, &
  dry_deposition, &
  chemistry_timestep
```

### Parameter Validation

```fortran
subroutine validate_config()
  if (nlat <= 0 .or. nlon <= 0) then
    call error_handler('Invalid grid dimensions')
  end if

  if (dx <= 0.0_r8 .or. dy <= 0.0_r8) then
    call error_handler('Invalid grid spacing')
  end if

  ! Additional validation...
end subroutine
```

## Extension Points

### Adding New Physics Modules

1. **Create module file**: `canopy_newphysics_mod.F90`
2. **Follow naming conventions**:
   - Module: `canopy_newphysics_mod`
   - Subroutines: `newphysics_init`, `newphysics_calc`, `newphysics_finalize`
3. **Update dependencies**: Add to Makefile and main program
4. **Add configuration**: Extend namelists as needed
5. **Include tests**: Create unit and integration tests

### Adding New I/O Formats

1. **Create I/O module**: `canopy_newformat_io_mod.F90`
2. **Implement interface**: Follow existing I/O patterns
3. **Update format options**: Extend `infmt_opt`/`outfmt_opt`
4. **Add validation**: Include format-specific checks

### Adding New Chemical Species

1. **Extend species list**: Update `canopy_canvars_mod`
2. **Add emission factors**: Update `canopy_bioparm_mod`
3. **Implement reactions**: Extend `canopy_chem_mod`
4. **Update I/O**: Include in output routines

## Documentation Standards

### Code Documentation

- **Module headers**: Purpose, dependencies, author, date
- **Subroutine headers**: Purpose, arguments, algorithm
- **Inline comments**: Complex algorithms and non-obvious code
- **Variable descriptions**: Units and physical meaning

### API Documentation

Generated automatically from Doxygen comments:

```fortran
!> @brief Calculate leaf temperature using energy balance
!> @details This subroutine solves the leaf energy balance equation
!>          to determine leaf temperature for sunlit and shaded leaves
!> @param[in] ppfd Photosynthetic photon flux density (μmol/m²/s)
!> @param[in] tair Air temperature (K)
!> @param[out] tleaf Leaf temperature (K)
!> @author Model Development Team
!> @date 2024
subroutine calc_leaf_temperature(ppfd, tair, tleaf)
```

This architecture provides a solid foundation for model development, maintenance, and extension while ensuring computational efficiency and scientific accuracy.
