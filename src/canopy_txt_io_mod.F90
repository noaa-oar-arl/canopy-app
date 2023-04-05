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
        integer k, loc

        if (infmt_opt .eq. 1) then !only output text with 1D input

            write(*,*)  'Writing Text Output'
            write(*,*)  '-------------------------------'

            if (ifcanwind) then
                write(*,*)  'Writing canopy wind output'
                write(*,*)  '-------------------------------'
! ... save as text file
                open(10, file=TRIM(TXTPREFX)//'_output_canopy_wind.txt')
                write(10, '(a30, f6.1, a2)') 'Reference height, h: ', href_set, 'm'
                write(10, '(a30, i6)') 'Number of model layers: ', modlays
                write(10, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'WS (m/s)'
                do loc=1, nlat*nlon
                    do k=1, modlays
                        write(10, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                            zk(k), canWIND(loc, k)
                    end do
                end do
            end if

! ... save as text file
            if (ifcanwaf) then
                write(*,*)  'Writing canopy WAF output'
                write(*,*)  '-------------------------------'
                open(11, file=TRIM(TXTPREFX)//'_output_waf.txt')
                write(11, '(a30, f6.1)') 'Reference height, h: ', href_set, 'm'
                write(11, '(a8, a9, a19, a19, a11)') 'Lat', 'Lon', 'Canopy height (m)', 'Flame height (m)', 'WAF'
                do loc=1, nlat*nlon
                    write(11, '(f8.2, f9.2, f19.2, f19.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                        variables(loc)%fh, flameh(loc), waf(loc)
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
                    do k=1, modlays
                        write(12, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                            zk(k), Kz(loc,k)
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
                    do k=1, modlays
                        write(13, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                            zk(k), rjcf(loc,k)
                    end do
                end do
            end if
            if (ifcanbio) then
                write(*,*)  'Writing biogenic emissions'
                write(*,*)  '-------------------------------'
! ... save as text file
                open(13, file=TRIM(TXTPREFX)//'_output_bio.txt')
                write(13, '(a30, f6.1, a2)') 'Reference height, h: ', href_set, 'm'
                write(13, '(a30, i6)') 'Number of model layers: ', modlays
                write(13, '(a8, a9, a12, a15, a15, a15, a15, a15, a15, a15, a15, a15, a15,             &
                & a15, a15, a15, a15, a15, a15, a15, a15, a15)') 'Lat', 'Lon', 'Height (m)',  &
                    'emi_isop', 'emi_myrc', 'emi_sabi', 'emi_limo', 'emi_care', 'emi_ocim',     &
                    'emi_bpin', 'emi_apin', 'emi_mono', 'emi_farn', 'emi_cary', 'emi_sesq',     &
                    'emi_mbol', 'emi_meth', 'emi_acet', 'emi_co',   'emi_bvoc', 'emi_svoc', 'emi_ovoc'
                do loc=1, nlat*nlon
                    do k=1, modlays
                        write(13, '(f8.2, f9.2, f12.2, es15.7, es15.7, es15.7, es15.7, es15.7, es15.7, &
                        &        es15.7, es15.7, es15.7, es15.7, es15.7, es15.7, es15.7, es15.7,    &
                        &        es15.7, es15.7, es15.7, es15.7, es15.7)')  &
                            variables(loc)%lat, variables(loc)%lon, &
                            zk(k), emi_isop(loc,k), emi_myrc(loc,k), emi_sabi(loc,k), emi_limo(loc,k), &
                            emi_care(loc,k), emi_ocim(loc,k), emi_bpin(loc,k), emi_apin(loc,k),        &
                            emi_mono(loc,k), emi_farn(loc,k), emi_cary(loc,k), emi_sesq(loc,k),        &
                            emi_mbol(loc,k), emi_meth(loc,k), emi_acet(loc,k), emi_co(loc,k),          &
                            emi_bvoc(loc,k), emi_svoc(loc,k), emi_ovoc(loc,k)
                    end do
                end do
            end if
        end if

    END SUBROUTINE write_txt

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

END MODULE canopy_txt_io_mod
