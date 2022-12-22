SUBROUTINE canopy_write_ncf (OUTPREFX)

!-------------------------------------------------------------------------------
! Name:     Output netCDF File
! Purpose:  Create a netCDF file of Canopy-App Output
! Revised:  23 Dec 2022  Original version.  (P.C Campbell Spero)
!-------------------------------------------------------------------------------

    USE canopy_canopts_mod
    USE canopy_canmet_mod
    USE canopy_coord_mod
    USE canopy_canvars_mod
    USE netcdf

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT( IN )  :: OUTPREFX
    CHARACTER(LEN=256)              :: fl
    CHARACTER(LEN=16),  PARAMETER   :: pname      = 'CANOPY_WRITE_NCF'
    INTEGER                         :: cdfid_m, rcode
    INTEGER                         :: dim_nx
    INTEGER                         :: dim_ny
    INTEGER                         :: dim_nz
    INTEGER                         :: dim_time
    INTEGER                         :: dim_timestr
    INTEGER,  SAVE,     ALLOCATABLE :: id_fld     ( : )
    INTEGER,  SAVE                  :: it         = 0
    INTEGER,            PARAMETER   :: len_time   = 1
    INTEGER                         :: n
    INTEGER                         :: nn
    INTEGER                         :: ntot
    INTEGER                         :: nvars
    CHARACTER(LEN=32)               :: var

!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

    CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CREATING DIMENSION FOR ', a, &
    & /, 1x, '***   IN FILE ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9200 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR DEFINING VARIABLE ', a, &
    & /, 1x, '***   IN FILE ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CREATING ATTRIBUTE FOR', a, &
    & /, 1x, '***   IN FILE ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9350 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR ENDING DEFINITIONS ', &
    & /, 1x, '***   IN FILE ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR WRITING VARIABLE ', a, &
    & /, 1x, '***   TO FILE ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9500 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CREATING NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9700 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CLOSING NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"


    if (infmt_opt .eq. 0) then !Input format is 2D and output is written to 2D

!-------------------------------------------------------------------------------
! Allocate necessary variables.
!-------------------------------------------------------------------------------

        it = it + 1

        nvars = nfld2dxy + nfld2dxyt + nfld3dxyzt

        IF ( .NOT. ALLOCATED ( id_fld ) ) ALLOCATE ( id_fld ( nvars ) )

!-----------------------------------------------------------------------------
! Create netCDF file.
!-----------------------------------------------------------------------------

        fl = TRIM(OUTPREFX)//trim('.nc')

        rcode = nf90_create (fl, nf90_clobber, cdfid_m)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9500) TRIM(pname), TRIM(fl), TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

!-----------------------------------------------------------------------------
! Set up dimensions.
!-----------------------------------------------------------------------------

        var = "time"
        rcode = nf90_def_dim (cdfid_m, TRIM(var), nf90_unlimited, dim_time)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        var = "timestr"
        rcode = nf90_def_dim (cdfid_m, TRIM(var), len_time, dim_timestr)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        var = "west_east"
        rcode = nf90_def_dim (cdfid_m, TRIM(var), nlon, dim_nx)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        var = "south_north"
        rcode = nf90_def_dim (cdfid_m, TRIM(var), nlat, dim_ny)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        var = "levels"
        rcode = nf90_def_dim (cdfid_m, TRIM(var), modlays, dim_nz)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        !-----------------------------------------------------------------------------
        ! Define variables that will populate the file.
        !-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Time-independent 2d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld2dxy
            var = TRIM(fld2dxy(n)%fldname)
            rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                (/ dim_nx, dim_ny /), id_fld(n))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
        ENDDO
        ntot = nfld2dxy

!-------------------------------------------------------------------------------
! Time-varying 2d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld2dxyt
            nn = ntot + n
            var = TRIM(fld2dxyt(n)%fldname)
            rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                (/ dim_nx, dim_ny, dim_time /), id_fld(nn))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
        ENDDO
        ntot = ntot + nfld2dxyt

!-------------------------------------------------------------------------------
! Time-varying 3d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld3dxyzt
            nn = ntot + n
            var = TRIM(fld3dxyzt(n)%fldname)
            rcode = nf90_def_var (cdfid_m, TRIM(var), nf90_float,  &
                (/ dim_nx, dim_ny, dim_nz, dim_time /), id_fld(nn))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9200) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
        ENDDO
        ntot = ntot + nfld3dxyzt

        !-----------------------------------------------------------------------------
        ! Define global attributes.
        !-----------------------------------------------------------------------------

        CALL canopy_outncfglobal (cdfid_m, fl)

        !-----------------------------------------------------------------------------
        ! Define attributes for the variables.
        !-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Time-independent 2d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld2dxy
            var = TRIM(fld2dxy(n)%fldname)
            rcode = nf90_put_att (cdfid_m, id_fld(n), 'long_name',  &
                TRIM(fld2dxy(n)%long_name))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
            rcode = nf90_put_att (cdfid_m, id_fld(n), 'units', TRIM(fld2dxy(n)%units))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
        ENDDO
        ntot = nfld2dxy

!-------------------------------------------------------------------------------
! Time-varying 2d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld2dxyt
            nn = ntot + n
            var = TRIM(fld2dxyt(n)%fldname)
            rcode = nf90_put_att (cdfid_m, id_fld(nn), 'long_name',  &
                TRIM(fld2dxyt(n)%long_name))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
            rcode = nf90_put_att (cdfid_m, id_fld(nn), 'units',  &
                TRIM(fld2dxyt(n)%units))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
        ENDDO
        ntot = ntot + nfld2dxyt

!-------------------------------------------------------------------------------
! Time-varying 3d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld3dxyzt
            nn = ntot + n
            var = TRIM(fld3dxyzt(n)%fldname)
            rcode = nf90_put_att (cdfid_m, id_fld(nn), 'long_name',  &
                TRIM(fld3dxyzt(n)%long_name))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
            rcode = nf90_put_att (cdfid_m, id_fld(nn), 'units',  &
                TRIM(fld3dxyzt(n)%units))
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
        ENDDO
        ntot = ntot + nfld3dxyzt

!-------------------------------------------------------------------------------
! Take file out of "define mode" so that variables may be written to it.
!-------------------------------------------------------------------------------

        rcode = nf90_enddef (cdfid_m)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9350) TRIM(pname), TRIM(fl), TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

!-------------------------------------------------------------------------------
! Assign pointer variables to the respective arrays from canopy-app.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Time-independent 2d fields at cell centers.
!-------------------------------------------------------------------------------

        g_lat%fld = variables_2d%lat
        g_lon%fld = variables_2d%lon

!-------------------------------------------------------------------------------
! Time-varying 2d fields at cell centers.
!-------------------------------------------------------------------------------

        c_waf%fld = waf_2d

!-------------------------------------------------------------------------------
! Time-varying 3d fields at cell centers.
!-------------------------------------------------------------------------------
        c_canwind%fld = canWIND_3d
        c_Kz%fld      = Kz_3d
        c_rjcf%fld    = rjcf_3d

!-------------------------------------------------------------------------------
! Write variables.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Time-independent 2d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld2dxy
            var = TRIM(fld2dxy(n)%fldname)
            rcode = nf90_put_var (cdfid_m, id_fld(n), fld2dxy(n)%fld)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
        ENDDO

        ntot = nfld2dxy

!-------------------------------------------------------------------------------
! Time-varying 2d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld2dxyt
            nn = ntot + n
            var = TRIM(fld2dxyt(n)%fldname)
            rcode = nf90_put_var (cdfid_m, id_fld(nn), fld2dxyt(n)%fld,  &
                start = (/ 1, 1, it /),  &
                count = (/ fld2dxyt(n)%iend(1), &
                fld2dxyt(n)%iend(2), 1 /) )
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
        ENDDO
        ntot = ntot + nfld2dxyt

!-------------------------------------------------------------------------------
! Time-independent 3d fields at cell centers.
!-------------------------------------------------------------------------------

        DO n = 1, nfld3dxyzt
            nn = ntot + n
            var = TRIM(fld3dxyzt(n)%fldname)
            rcode = nf90_put_var (cdfid_m, id_fld(nn), fld3dxyzt(n)%fld,  &
                start = (/ 1, 1, 1, it /),  &
                count = (/ fld3dxyzt(n)%iend(1), &
                fld3dxyzt(n)%iend(2), &
                fld3dxyzt(n)%iend(3), 1 /) )
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9400) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit (2)
            ENDIF
        ENDDO
        ntot = ntot + nfld3dxyzt

!-------------------------------------------------------------------------------
! Gracefully close output files.
!-------------------------------------------------------------------------------

        rcode = nf90_close (cdfid_m)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9700) TRIM(pname), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

    end if !infmt_opt

END SUBROUTINE canopy_write_ncf
