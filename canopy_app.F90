program canopy_app
!
!  The main application
!
!  History:
!    Prototype: Patrick C. Campbell, 06/2022
!    Revised  : PCC (10/2022)
!-------------------------------------------------------------
    use canopy_files_mod   !main canopy input files

    implicit none

!-------------------------------------------------------------------------------
! Read user options from namelist.
!-------------------------------------------------------------------------------

    call canopy_readnml

!-------------------------------------------------------------------------------
! Allocate necessary variables.
!-------------------------------------------------------------------------------

    call canopy_alloc

!-------------------------------------------------------------------------------
! Initialize canopy variables
!-------------------------------------------------------------------------------

    call canopy_init

!-------------------------------------------------------------------------------
! Read met/sfc gridded model input file (currently 1D TXT or 1D/2D NETCDF).
!-------------------------------------------------------------------------------

    call canopy_check_input(file_vars(1))

!-------------------------------------------------------------------------------
! Main canopy model calculations.
!-------------------------------------------------------------------------------

    call canopy_calcs

!-------------------------------------------------------------------------------
! Allocate and initialize 2D/3D NetCDF output data structures
!-------------------------------------------------------------------------------
!TODO--->
!    call canopy_outncf_alloc
!    call canopy_outncf_init

!-------------------------------------------------------------------------------
! Write model output of canopy model calculations.
!-------------------------------------------------------------------------------

    call canopy_write_txt(file_out(1))

    !TODO--->
!    if (infmt_opt .eq. 0) then !Input format is 2D, then add and 2D NetCDF output
!
!        call canout_write_ncf(file_out(1))
!
!    end if

!-------------------------------------------------------------------------------
! Dellocate necessary variables.
!-------------------------------------------------------------------------------

    call canopy_dealloc

end program canopy_app
