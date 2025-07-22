!> \file canopy_files_mod.F90
!! \brief File Management Module
!! \details This module contains FORTRAN units and file names for the canopy model.
!! It manages input and output file specifications including variable files,
!! canopy profile files, output files, and the namelist file.
!!
!! \author Patrick C. Campbell
!! \date July 2022
!!
!! \version
!! - 15 Jul 2022: Original version (P. C. Campbell)
!! - 30 Nov 2023: Added supplementary canopy profile, file_canvars (P.C. Campbell)

!> \defgroup files_mod File Management Module
!! \brief Module for file path and unit management
!! \{

MODULE canopy_files_mod

    IMPLICIT NONE

!> \defgroup file_units File Units and Limits
!! \brief File unit numbers and array size limits
!! \{

    !> \brief NetCDF file identifier for main input files
    INTEGER                       :: cdfid_m
    !> \brief Maximum number of input files allowed
    INTEGER,            PARAMETER :: max_mm       = 10000
    !> \brief FORTRAN unit number for namelist file
    INTEGER,            PARAMETER :: iutnml       =  8

!> \}

!> \defgroup file_arrays File Path Arrays
!! \brief Arrays storing file paths for different input/output types
!! \{

    !> \brief Array of main variable input file paths
    CHARACTER(LEN=256)            :: file_vars    ( max_mm )
    !> \brief Array of canopy variable input file paths
    CHARACTER(LEN=256)            :: file_canvars ( max_mm )
    !> \brief Array of output file paths
    CHARACTER(LEN=256)            :: file_out     ( 1 )

!> \}

!> \defgroup file_paths Fixed File Paths
!! \brief Fixed file paths for configuration files
!! \{

    !> \brief Path to the namelist configuration file
    CHARACTER(LEN=*), PARAMETER   :: file_nml     = 'input/namelist.canopy'

!> \}

!> \}

END MODULE canopy_files_mod
