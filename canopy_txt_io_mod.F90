MODULE canopy_txt_io_mod

!-------------------------------------------------------------------------------
! Name:     TXT IO
! Purpose:  Contains routines to read met/sfc model text output.
! Revised:  03 Oct 2022  Original version  (P.C. Campbell)
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE read_txt(TXTFILE)

        USE canopy_coord_mod
        USE canopy_canmet_mod

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT( IN )  :: TXTFILE

        !Local variables
        integer i0, loc

        ! ... read met/sfc input variables from text file
        open(8,  file=TXTFILE,  status='old')
        i0 = 0
        read(8,*,iostat=i0)  ! skip headline
        do loc=1, nlat*nlon
            read(8, *) variables(loc)
        end do
        close(8)

    END SUBROUTINE read_txt

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


    SUBROUTINE write_txt(TXTPREFX)

        USE canopy_coord_mod
        USE canopy_canopts_mod
        USE canopy_canmet_mod
        USE canopy_canvars_mod

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT( IN )  :: TXTPREFX

        !Local variables
        integer i, loc

        if (ifcanwind) then
            write(*,*)  'Writing canopy wind output'
            write(*,*)  '-------------------------------'
! ... save as text file
            open(10, file=TRIM(TXTPREFX)//'_output_canopy_wind.txt')
            write(10, '(a30, f6.1, a2)') 'Reference height, h: ', href_set, 'm'
            write(10, '(a30, i6)') 'Number of model layers: ', modlays
            write(10, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'WS (m/s)'
            do loc=1, nlat*nlon
                do i=1, modlays
                    write(10, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                        zk(i), canWIND(i, loc)
                end do
            end do
        end if

! ... save as text file
        if (ifcanwaf) then
            write(*,*)  'Writing canopy WAF output'
            write(*,*)  '-------------------------------'
            open(11, file=TRIM(TXTPREFX)//'_output_waf.txt')
            write(11, '(a30, f6.1)') 'Reference height, h: ', href_set, 'm'
            write(11, '(a8, a9, a19, a11)') 'Lat', 'Lon', 'Canopy height (m)', 'WAF'
            do loc=1, nlat*nlon
                write(11, '(f8.2, f9.2, f19.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, hcmref, waf(loc)
            end do
        end if

        if (ifcaneddy) then
            write(*,*)  'Writing canopy eddy diffusivity scaling values'
            write(*,*)  '-------------------------------'
! ... save as text file
            open(12, file=TRIM(TXTPREFX)//'_output_eddy_Kz.txt')
            write(12, '(a30, f6.1, a2)') 'Reference height, h: ', href_set, 'm'
            write(12, '(a30, i6)') 'Number of model layers: ', modlays
            write(12, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'Kz'
            do loc=1, nlat*nlon
                do i=1, modlays
                    write(12, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                        zk(i), Kz(i, loc)
                end do
            end do
        end if

        if (ifcanphot) then
            write(*,*)  'Writing canopy photolysis correction factors'
            write(*,*)  '-------------------------------'
! ... save as text file
            open(13, file=TRIM(TXTPREFX)//'_output_phot.txt')
            write(13, '(a30, f6.1, a2)') 'Reference height, h: ', href_set, 'm'
            write(13, '(a30, i6)') 'Number of model layers: ', modlays
            write(13, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'rjcf'
            do loc=1, nlat*nlon
                do i=1, modlays
                    write(13, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                        zk(i), rjcf(i, loc)
                end do
            end do
        end if

    END SUBROUTINE write_txt

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

END MODULE canopy_txt_io_mod
