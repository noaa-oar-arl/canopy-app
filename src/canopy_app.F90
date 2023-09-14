program canopy_app
!
!  The main application
!
!  History:
!    Prototype: Patrick C. Campbell, 06/2022
!    Revised  : PCC (10/2022)
!    21 Aug 2023 Adding multiple timesteps (P.C. Campbell)
!-------------------------------------------------------------
    use canopy_date_mod   ! main canopy date module
    use canopy_files_mod  ! main canopy input files
    use canopy_coord_mod  ! main canopy coordinates
#ifdef NETCDF
    use canopy_ncf_io_mod
#endif

    implicit none

    CHARACTER(LEN=24)                 :: time_next  ! YYYY-MM-DD-HH:MM:SS.SSSS
    CHARACTER(LEN=24)                 :: time_now   ! YYYY-MM-DD-HH:MM:SS.SSSS
    integer   :: nn
    CHARACTER(LEN=24)                 :: nn_string
!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

    CHARACTER(LEN=256), PARAMETER :: f100 = "(//, 1x, 78('~'), &
    & /,  1x, '~~~ Processing canopy-app for time = ', a,  &
    & /,  1x, 78('~'), /)"

    !   CHARACTER(LEN=256), PARAMETER :: f200 = "(//, 1x, 78('~'), &
    !   & /,  1x, '~~~ Metadata summary', &
    !   & /,  1x, 78('~'), /)"

!-------------------------------------------------------------------------------


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
! Allocate and initialize 2D/3D NetCDF output data structures
!-------------------------------------------------------------------------------

#ifdef NETCDF
    call canopy_outncf_alloc

    call canopy_outncf_init
#endif

!-------------------------------------------------------------------------------
! Loop over time to get input, calculate canopy fields, and write output.
!-------------------------------------------------------------------------------
    time_now = time_start
    if(ntime.le.0) ntime=999999999 ! assign a large number
    timeloop: DO nn=1,ntime

        WRITE (*,f100) time_now
!-------------------------------------------------------------------------------
! Read met/sfc gridded model input file (currently 1D TXT or 1D/2D NETCDF).
!-------------------------------------------------------------------------------

#ifdef NETCDF
        call canopy_check_input(file_vars(nn))
#else
        call canopy_read_txt(file_vars(nn))
#endif

!-------------------------------------------------------------------------------
! Main canopy model calculations.
!-------------------------------------------------------------------------------

        call canopy_calcs

!-------------------------------------------------------------------------------
! Write model output of canopy model calculations.
!-------------------------------------------------------------------------------
        write(nn_string,*) nn-1
!        call canopy_write_txt(trim(file_out(1)) // '_' // trim(time_now))
        if (nn.lt.10) then
            call canopy_write_txt((trim(file_out(1)) // '_t00' // ADJUSTL(nn_string)), &
                time_now)
        else if (nn.ge.10.and.nn.lt.100) then
            call canopy_write_txt((trim(file_out(1)) // '_t0' // ADJUSTL(nn_string)), &
                time_now)
        else
            call canopy_write_txt((trim(file_out(1)) // '_t' // ADJUSTL(nn_string)), &
                time_now)
        end if

#ifdef NETCDF
        !only output if 2D input NCF is used
        call canopy_write_ncf(trim(file_out(1)))
#endif

! Update new date and time

        CALL geth_newdate (time_next, time_now, time_intvl)
        time_now = time_next

    ENDDO timeloop

!-------------------------------------------------------------------------------
! Dellocate necessary variables.
!-------------------------------------------------------------------------------

    call canopy_dealloc

!-------------------------------------------------------------------------------
! Close output NetCDF files (if necessary).
!-------------------------------------------------------------------------------

!#ifdef NETCDF
!        call canopy_close_files(file_out(1))
!#endif

    WRITE (*,'(//, a)') 'Canopy-App Finished Normally'

end program canopy_app
