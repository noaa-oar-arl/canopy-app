!> \file canopy_coord_mod.F90
!! \brief Canopy Coordinate and Domain Module
!! \details This module contains coordinate and domain descriptions for the canopy model
!! including time and geographic domain specifications that are read from files and/or
!! user namelist input.
!!
!! \author Patrick C. Campbell
!! \date October 2022

!> \defgroup coord_mod Canopy Coordinate and Domain Module
!! \brief Module for time and spatial domain specifications
!! \{

MODULE canopy_coord_mod
    use canopy_const_mod, ONLY: rk
    IMPLICIT NONE

!> \defgroup time_domain Time Domain Variables
!! \brief Variables defining the temporal domain of the simulation
!! \{

    !> \brief Simulation start time in YYYY-MM-DD-HH:MM:SS.SSSS format
    CHARACTER(LEN=24)  :: time_start
    !> \brief Simulation end time in YYYY-MM-DD-HH:MM:SS.SSSS format
    CHARACTER(LEN=24)  :: time_end
    !> \brief Time interval for input/output [seconds]
    INTEGER            :: time_intvl
    !> \brief Number of model timesteps
    integer            :: ntime

!> \}

!> \defgroup spatial_domain Spatial Domain Variables
!! \brief Variables defining the spatial domain of the simulation
!! \{

    !> \brief Length of latitude coordinate (number of latitude points)
    integer            :: nlat
    !> \brief Length of longitude coordinate (number of longitude points)
    integer            :: nlon

!> \}

!> \defgroup vertical_domain Vertical Domain Variables
!! \brief Variables defining the vertical domain of the canopy model
!! \{

    !> \brief Number of total above and below canopy model layers
    integer            :: modlays
    !> \brief Model above and below canopy vertical resolution [m]
    real(rk)           :: modres

!> \}

!> \}

END MODULE canopy_coord_mod
