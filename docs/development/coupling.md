# Model Coupling

This document describes how to couple the Canopy-App model with other atmospheric, land surface, and Earth system models. It provides guidance for developers implementing model coupling interfaces.

## Coupling Overview

### Types of Coupling

The Canopy-App model supports several coupling approaches:

1. **One-way coupling**: Canopy-App receives inputs from another model
2. **Two-way coupling**: Bidirectional exchange of information
3. **Component coupling**: Integration as a module within a larger system
4. **Framework coupling**: Using coupling frameworks like ESMF or OASIS

### Coupling Variables

#### Variables Received from Atmosphere Models

| Variable | Units | Description | Frequency |
|----------|-------|-------------|-----------|
| `temp_air` | K | Air temperature | Every timestep |
| `qv_air` | kg/kg | Specific humidity | Every timestep |
| `press_air` | Pa | Air pressure | Every timestep |
| `wind_u` | m/s | U-component wind | Every timestep |
| `wind_v` | m/s | V-component wind | Every timestep |
| `srad_down` | W/m² | Downward shortwave radiation | Every timestep |
| `lrad_down` | W/m² | Downward longwave radiation | Every timestep |
| `precip` | mm/s | Precipitation rate | Every timestep |

#### Variables Sent to Atmosphere Models

| Variable | Units | Description | Frequency |
|----------|-------|-------------|-----------|
| `sensible_heat` | W/m² | Sensible heat flux | Every timestep |
| `latent_heat` | W/m² | Latent heat flux | Every timestep |
| `momentum_flux_u` | N/m² | U-momentum flux | Every timestep |
| `momentum_flux_v` | N/m² | V-momentum flux | Every timestep |
| `co2_flux` | μmol/m²/s | CO₂ flux | Every timestep |
| `albedo` | - | Surface albedo | As needed |
| `roughness_length` | m | Aerodynamic roughness | As needed |
| `biogenic_emissions` | μg/m²/s | Biogenic VOC emissions | Every timestep |

## Coupling Interfaces

### Generic Coupling Interface

The model provides a generic coupling interface that can be adapted for different frameworks:

```fortran
module canopy_coupling_mod
  use canopy_const_mod, only: r8
  implicit none
  private

  ! Public coupling procedures
  public :: canopy_coupling_init
  public :: canopy_coupling_import
  public :: canopy_coupling_export
  public :: canopy_coupling_finalize

  ! Coupling state
  type :: coupling_state_type
    logical :: is_coupled
    integer :: coupling_timestep
    character(len=256) :: coupling_method
  end type

  type(coupling_state_type), save :: coupling_state

contains

  subroutine canopy_coupling_init(method, timestep)
    character(len=*), intent(in) :: method
    integer, intent(in) :: timestep

    coupling_state%coupling_method = method
    coupling_state%coupling_timestep = timestep
    coupling_state%is_coupled = .true.

    select case (trim(method))
    case ('ESMF')
      call init_esmf_coupling()
    case ('OASIS')
      call init_oasis_coupling()
    case ('MCT')
      call init_mct_coupling()
    case ('INTERNAL')
      call init_internal_coupling()
    case default
      call error_handler('Unknown coupling method: ' // trim(method))
    end select
  end subroutine

end module
```

### ESMF Coupling

Earth System Modeling Framework integration:

```fortran
module canopy_esmf_coupling_mod
  use ESMF
  use canopy_canvars_mod
  implicit none
  private

  public :: canopy_esmf_init, canopy_esmf_run, canopy_esmf_finalize

  ! ESMF components
  type(ESMF_GridComp) :: canopy_comp
  type(ESMF_State) :: import_state, export_state
  type(ESMF_Grid) :: grid

contains

  subroutine canopy_esmf_init(gridcomp, importstate, exportstate, clock, rc)
    type(ESMF_GridComp), intent(inout) :: gridcomp
    type(ESMF_State), intent(inout) :: importstate
    type(ESMF_State), intent(inout) :: exportstate
    type(ESMF_Clock), intent(in) :: clock
    integer, intent(out) :: rc

    ! Create grid
    grid = ESMF_GridCreateNoPeriDim(minIndex=[1,1], maxIndex=[nlon,nlat], &
                                   coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)

    ! Add coordinates
    call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    ! Create fields and add to states
    call create_import_fields(importstate, grid, rc)
    call create_export_fields(exportstate, grid, rc)

    ! Initialize canopy model
    call canopy_init()

    rc = ESMF_SUCCESS
  end subroutine

  subroutine canopy_esmf_run(gridcomp, importstate, exportstate, clock, rc)
    type(ESMF_GridComp), intent(inout) :: gridcomp
    type(ESMF_State), intent(inout) :: importstate
    type(ESMF_State), intent(inout) :: exportstate
    type(ESMF_Clock), intent(in) :: clock
    integer, intent(out) :: rc

    ! Get import data
    call get_import_data(importstate, rc)

    ! Run canopy model
    call canopy_calcs(1)  ! Assuming single timestep

    ! Set export data
    call set_export_data(exportstate, rc)

    rc = ESMF_SUCCESS
  end subroutine

end module
```

### OASIS Coupling

OASIS coupler integration:

```fortran
module canopy_oasis_coupling_mod
  use mod_oasis
  use canopy_canvars_mod
  implicit none
  private

  public :: init_oasis_coupling, oasis_exchange, finalize_oasis_coupling

  ! OASIS variables
  integer :: comp_id
  integer :: partition_id
  integer :: var_ids(20)  ! Assuming max 20 coupling variables
  integer :: n_vars

contains

  subroutine init_oasis_coupling()
    integer :: ierror, local_comm
    character(len=6) :: comp_name = 'canopy'

    ! Initialize OASIS
    call oasis_init_comp(comp_id, comp_name, ierror)

    ! Get local communicator
    call oasis_get_localcomm(local_comm, ierror)

    ! Define partition
    call define_partition()

    ! Define coupling variables
    call define_coupling_variables()

    ! End definition phase
    call oasis_enddef(ierror)
  end subroutine

  subroutine oasis_exchange(step)
    integer, intent(in) :: step
    integer :: ierror

    ! Receive data from atmosphere
    call oasis_get(var_ids(1), step, temp_air, ierror)  ! Temperature
    call oasis_get(var_ids(2), step, qv_air, ierror)    ! Humidity
    call oasis_get(var_ids(3), step, wind_u, ierror)    ! U-wind
    call oasis_get(var_ids(4), step, wind_v, ierror)    ! V-wind

    ! Run canopy model
    call canopy_calcs(step)

    ! Send data to atmosphere
    call oasis_put(var_ids(5), step, sensible_heat, ierror)
    call oasis_put(var_ids(6), step, latent_heat, ierror)
    call oasis_put(var_ids(7), step, momentum_flux_u, ierror)
    call oasis_put(var_ids(8), step, momentum_flux_v, ierror)
  end subroutine

end module
```

## Specific Model Couplings

### WRF Coupling

Weather Research and Forecasting model integration:

```fortran
module canopy_wrf_coupling_mod
  implicit none
  private

  public :: canopy_wrf_init, canopy_wrf_driver

contains

  subroutine canopy_wrf_init(ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte)
    integer, intent(in) :: ids, ide, jds, jde, kds, kde
    integer, intent(in) :: ims, ime, jms, jme, kms, kme
    integer, intent(in) :: its, ite, jts, jte, kts, kte

    ! Initialize canopy model with WRF grid dimensions
    call canopy_init_grid(its, ite, jts, jte, kts, kte)
    call canopy_init()
  end subroutine

  subroutine canopy_wrf_driver(u, v, th, qv, p, rho, dt, &
                              dx, dy, dz, z, &
                              tsk, hfx, qfx, lh, &
                              its, ite, jts, jte, kts, kte)
    ! WRF state variables
    real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: u, v, th, qv, p, rho
    real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: dz, z
    real, dimension(ims:ime,jms:jme), intent(in) :: tsk
    real, dimension(ims:ime,jms:jme), intent(inout) :: hfx, qfx, lh
    real, intent(in) :: dt, dx, dy
    integer, intent(in) :: its, ite, jts, jte, kts, kte

    ! Copy WRF data to canopy arrays
    call wrf_to_canopy_transfer(u, v, th, qv, p, rho, tsk, &
                               its, ite, jts, jte, kts, kte)

    ! Run canopy model
    call canopy_calcs(1)

    ! Copy canopy results back to WRF
    call canopy_to_wrf_transfer(hfx, qfx, lh, &
                               its, ite, jts, jte)
  end subroutine

end module
```

### CESM Coupling

Community Earth System Model integration:

```fortran
module canopy_cesm_coupling_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use seq_cdata_mod
  use seq_infodata_mod
  implicit none
  private

  public :: canopy_cesm_init, canopy_cesm_run, canopy_cesm_final

contains

  subroutine canopy_cesm_init(EClock, cdata, x2l, l2x, NLFilename)
    type(seq_cdata), intent(inout) :: cdata
    type(mct_aVect), intent(inout) :: x2l, l2x
    character(len=*), optional, intent(in) :: NLFilename

    ! Initialize canopy model
    call canopy_init()

    ! Set up attribute vectors for coupling
    call setup_coupling_vectors(x2l, l2x)
  end subroutine

  subroutine canopy_cesm_run(EClock, cdata, x2l, l2x)
    type(seq_cdata), intent(inout) :: cdata
    type(mct_aVect), intent(inout) :: x2l, l2x

    integer :: lsize

    ! Get data from coupler
    lsize = mct_aVect_lsize(x2l)
    call receive_atm_data(x2l, lsize)

    ! Run canopy model
    call canopy_calcs(1)

    ! Send data to coupler
    call send_lnd_data(l2x, lsize)
  end subroutine

end module
```

### FV3 Coupling

NOAA's Finite-Volume Cubed-Sphere Dynamical Core:

```fortran
module canopy_fv3_coupling_mod
  use canopy_canvars_mod
  implicit none
  private

  public :: canopy_fv3_init, canopy_fv3_run

contains

  subroutine canopy_fv3_init(nlon, nlat, nlev)
    integer, intent(in) :: nlon, nlat, nlev

    ! Initialize with FV3 grid
    call canopy_init_dimensions(nlon, nlat, nlev)
    call canopy_init()
  end subroutine

  subroutine canopy_fv3_run(temp, qv, u, v, ps, &
                           hflx, qflx, uflx, vflx, &
                           nlon, nlat, nlev, dt)
    integer, intent(in) :: nlon, nlat, nlev
    real(r8), intent(in) :: dt
    real(r8), dimension(nlon,nlat,nlev), intent(in) :: temp, qv, u, v
    real(r8), dimension(nlon,nlat), intent(in) :: ps
    real(r8), dimension(nlon,nlat), intent(out) :: hflx, qflx, uflx, vflx

    ! Transfer FV3 data to canopy
    call fv3_to_canopy(temp, qv, u, v, ps, nlon, nlat, nlev)

    ! Run canopy model
    call canopy_calcs(1)

    ! Transfer canopy results to FV3
    call canopy_to_fv3(hflx, qflx, uflx, vflx, nlon, nlat)
  end subroutine

end module
```

## Data Exchange Patterns

### Synchronous Coupling

Both models run with the same timestep:

```fortran
subroutine synchronous_coupling_step(atm_dt)
  real(r8), intent(in) :: atm_dt

  ! Receive atmospheric data
  call receive_from_atmosphere()

  ! Run canopy model for same timestep
  call canopy_calcs_dt(atm_dt)

  ! Send surface fluxes back
  call send_to_atmosphere()
end subroutine
```

### Asynchronous Coupling

Models run with different timesteps:

```fortran
subroutine asynchronous_coupling_step(atm_dt, cnp_dt)
  real(r8), intent(in) :: atm_dt, cnp_dt
  integer :: n_substeps

  ! Calculate number of canopy substeps
  n_substeps = nint(atm_dt / cnp_dt)

  ! Receive atmospheric data (valid for entire atm_dt)
  call receive_from_atmosphere()

  ! Run multiple canopy timesteps
  do istep = 1, n_substeps
    call canopy_calcs_dt(cnp_dt)

    ! Accumulate fluxes
    call accumulate_fluxes(istep, n_substeps)
  end do

  ! Send time-averaged fluxes
  call send_averaged_fluxes()
end subroutine
```

### Conservative Interpolation

For grid mismatches between models:

```fortran
module conservative_interpolation_mod
  implicit none
  private

  public :: setup_interpolation, interpolate_conservative

  type :: interp_weight_type
    integer :: n_src, n_dst
    integer, allocatable :: src_indices(:)
    integer, allocatable :: dst_indices(:)
    real(r8), allocatable :: weights(:)
  end type

contains

  subroutine setup_interpolation(src_grid, dst_grid, weights)
    type(grid_type), intent(in) :: src_grid, dst_grid
    type(interp_weight_type), intent(out) :: weights

    ! Calculate conservative interpolation weights
    ! This is a simplified version - real implementation would be more complex
    call calculate_overlap_areas(src_grid, dst_grid, weights)
  end subroutine

  subroutine interpolate_conservative(src_data, dst_data, weights)
    real(r8), intent(in) :: src_data(:)
    real(r8), intent(out) :: dst_data(:)
    type(interp_weight_type), intent(in) :: weights

    integer :: i

    ! Initialize destination
    dst_data = 0.0_r8

    ! Apply weights
    do i = 1, weights%n_weights
      dst_data(weights%dst_indices(i)) = dst_data(weights%dst_indices(i)) + &
        weights%weights(i) * src_data(weights%src_indices(i))
    end do
  end subroutine

end module
```

## Testing Coupled Systems

### Unit Testing Coupling Interfaces

```fortran
module test_coupling_mod
  use canopy_coupling_mod
  use canopy_test_framework
  implicit none

contains

  subroutine test_coupling_init()
    integer :: status

    ! Test coupling initialization
    call canopy_coupling_init('INTERNAL', 3600, status)
    call assert_equal(status, 0, 'Coupling initialization failed')

    ! Test coupling state
    call assert_true(coupling_state%is_coupled, 'Coupling not active')
    call assert_equal(coupling_state%coupling_timestep, 3600, 'Wrong timestep')
  end subroutine

  subroutine test_data_exchange()
    real(r8) :: temp_in(10,10) = 298.15_r8
    real(r8) :: flux_out(10,10)

    ! Import test data
    call canopy_coupling_import('temperature', temp_in)

    ! Run model
    call canopy_calcs(1)

    ! Export test data
    call canopy_coupling_export('sensible_heat', flux_out)

    ! Verify reasonable flux values
    call assert_range(flux_out, -200.0_r8, 800.0_r8, 'Unrealistic heat flux')
  end subroutine

end module
```

### Integration Testing

```bash
#!/bin/bash
# test_coupled_system.sh

# Test WRF coupling
echo "Testing WRF coupling..."
cd test/wrf_coupling
make clean && make
mpirun -np 4 ./wrf.exe
if [ $? -eq 0 ]; then
  echo "WRF coupling test PASSED"
else
  echo "WRF coupling test FAILED"
  exit 1
fi

# Test CESM coupling
echo "Testing CESM coupling..."
cd ../cesm_coupling
make clean && make
mpirun -np 8 ./cesm.exe
if [ $? -eq 0 ]; then
  echo "CESM coupling test PASSED"
else
  echo "CESM coupling test FAILED"
  exit 1
fi
```

## Best Practices

### Performance Considerations

1. **Minimize data copying**: Use pointers and references where possible
2. **Efficient interpolation**: Pre-compute interpolation weights
3. **Communication optimization**: Aggregate small messages
4. **Load balancing**: Ensure even distribution of computational work

### Error Handling

```fortran
subroutine robust_coupling_exchange(status, error_msg)
  integer, intent(out) :: status
  character(len=*), intent(out) :: error_msg

  status = 0
  error_msg = ''

  ! Receive data with error checking
  call receive_coupling_data(status)
  if (status /= 0) then
    error_msg = 'Failed to receive coupling data'
    return
  end if

  ! Validate received data
  call validate_coupling_data(status)
  if (status /= 0) then
    error_msg = 'Invalid coupling data received'
    return
  end if

  ! Run model with error handling
  call canopy_calcs_safe(status)
  if (status /= 0) then
    error_msg = 'Canopy model calculation failed'
    return
  end if

  ! Send results
  call send_coupling_data(status)
  if (status /= 0) then
    error_msg = 'Failed to send coupling data'
    return
  end if

end subroutine
```

### Coupling Validation

```fortran
subroutine validate_coupling_conservation()
  real(r8) :: mass_in, mass_out, energy_in, energy_out
  real(r8) :: mass_error, energy_error

  ! Calculate mass balance
  mass_in = sum(precip_rate) * dt
  mass_out = sum(evap_rate) * dt
  mass_error = abs(mass_in - mass_out) / max(mass_in, 1.0e-10_r8)

  ! Calculate energy balance
  energy_in = sum(net_radiation)
  energy_out = sum(sensible_heat + latent_heat)
  energy_error = abs(energy_in - energy_out) / max(energy_in, 1.0e-10_r8)

  ! Check conservation
  if (mass_error > 0.01_r8) then  ! 1% tolerance
    call warning_handler('Mass conservation error: ' // &
                         real_to_string(mass_error*100) // '%')
  end if

  if (energy_error > 0.05_r8) then  ! 5% tolerance
    call warning_handler('Energy conservation error: ' // &
                         real_to_string(energy_error*100) // '%')
  end if

end subroutine
```

This coupling framework provides a foundation for integrating the Canopy-App model with various Earth system models while maintaining computational efficiency and scientific accuracy.
