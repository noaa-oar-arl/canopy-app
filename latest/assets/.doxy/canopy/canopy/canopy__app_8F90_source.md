

# File canopy\_app.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_app.F90**](canopy__app_8F90.md)

[Go to the documentation of this file](canopy__app_8F90.md)


```Fortran



program canopy_app
    use canopy_date_mod    ! main canopy date module
    use canopy_files_mod   ! main canopy input files
    use canopy_coord_mod   ! main canopy coordinates

#ifdef NETCDF
    use canopy_ncf_io_mod
#endif

    implicit none

    CHARACTER(LEN=24)                 :: time_next
    CHARACTER(LEN=24)                 :: time_now
    CHARACTER(LEN=24)                 :: time_now_file

    integer   :: n,nn
!    CHARACTER(LEN=24)                 :: nn_string

    CHARACTER(LEN=256), PARAMETER :: f100 = "(//, 1x, 78('~'), &
    & /,  1x, '~~~ Processing canopy-app for time = ', a,  &
    & /,  1x, 78('~'), /)"




    call canopy_readnml

    call canopy_alloc

    call canopy_init

#ifdef NETCDF

    call canopy_outncf_alloc

    call canopy_outncf_init
#endif

    time_now = time_start
    if(ntime.le.0) ntime=999999999 ! assign a large number
    timeloop: DO nn=1,ntime

        WRITE (*,f100) time_now
#ifdef NETCDF

        call canopy_check_input(file_vars(nn),file_canvars(nn))
#else

        call canopy_read_txt(file_vars(nn),file_canvars(nn))
#endif

        call canopy_calcs(nn)

        time_now_file=time_now
        do n = 1, len(time_now_file)
            if (time_now_file(n:n) == ':') then
                time_now_file(n:n) = '-'
            end if
        end do

        call canopy_write_txt((trim(file_out(1)) // '_' // trim(time_now_file)),time_now)

!        if (nn.lt.10) then
!            call canopy_write_txt((trim(file_out(1)) // '_t00' // ADJUSTL(nn_string)), &
!                time_now)
!        else if (nn.ge.10.and.nn.lt.100) then
!            call canopy_write_txt((trim(file_out(1)) // '_t0' // ADJUSTL(nn_string)), &
!                time_now)
!        else
!            call canopy_write_txt((trim(file_out(1)) // '_t' // ADJUSTL(nn_string)), &
!                time_now)
!        end if

#ifdef NETCDF

        call canopy_write_ncf(trim(file_out(1)) // '_' // adjustl(time_now_file))
#endif

        CALL geth_newdate (time_next, time_now, time_intvl)
        time_now = time_next

    ENDDO timeloop

    call canopy_dealloc



    WRITE (*,'(//, a)') 'Canopy-App Finished Normally'

end program canopy_app
```


