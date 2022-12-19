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

    !Local variables
    integer ppos

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

!    TODO :: Add a separate intialize routine for gridded ncdf canopy variables
!            to fit ncf file, e.g.,, call canout_ncf_init

!-------------------------------------------------------------------------------
! Read met/sfc gridded model input file (currently text or 1D or 2D ncf).
!-------------------------------------------------------------------------------
! ... TODO: Read from txt or 1D/2D netcdf file
    ppos = scan(trim(file_vars(1)),".", BACK= .true.)
    if (trim(file_vars(1)(ppos:)).eq.".txt") then !TXT File
        call canopy_read_txt(file_vars(1))
    else if (trim(file_vars(1)(ppos:)).eq.".nc") then !NetCDF File
        call canopy_read_ncf(file_vars(1))
    else if (trim(file_vars(1)(ppos:)).eq.".ncf") then
        call canopy_read_ncf(file_vars(1))
    else if (trim(file_vars(1)(ppos:)).eq.".nc4") then
        call canopy_read_ncf(file_vars(1))
    else
        write(*,*)  'Error the file input type ',trim(file_vars(1)(ppos:)), &
            ' is not supported...exiting'
        call exit(2)
    end if   !File Input types

! ... TODO:  Set nlat = nlat_user (if txt) or nlat = nlat_file (if ncf)
! ... TODO:  Set nlon = nlat_user (if txt) or nlon = nlon_file (if ncf)

!-------------------------------------------------------------------------------
! Main canopy model calculations.
!-------------------------------------------------------------------------------

    call canopy_calcs

!-------------------------------------------------------------------------------
! Write model output of canopy model calculations.
!-------------------------------------------------------------------------------

! ... TODO:  NL condition for data write to txt or netcdf 1D or 2D

    call canopy_write_txt(file_out(1))

!-------------------------------------------------------------------------------
! Dellocate necessary variables.
!-------------------------------------------------------------------------------

    call canopy_dealloc

end program canopy_app
