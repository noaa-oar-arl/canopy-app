MODULE canopy_ncf_io_mod

    public canopy_outncf_alloc, canopy_outncf_init, canopy_outncfglobal, &
        canopy_read_ncf, canopy_write_ncf
    private

!-------------------------------------------------------------------------------
! Name:     NetCDF IO
! Purpose:  Contains routines to read NetCDF output.
! Revised:  31 Aug 2004  Original version provided by U. Houston.  (S.-B. Kim)
!           21 Jul 2005  Updated to include error-handling on NetCDF functions
!                        and to conform to NetCDF standard routines for
!                        retrieving data (to eliminate type mismatches).
!                        Updated formatting.  (T. Otte)
!           02 Aug 2005  Changed order of variable declarations in some
!                        subroutines to avoid compile failure on some
!                        machines.  (T. Otte)
!           20 Jun 2006  Removed unused variables.  Changed local variables
!                        FILE to FILENAME and DATA to DATAOUT to avoid
!                        conflicts with F90 keywords.  (T. Otte)
!           19 Apr 2007  Added new routine GET_VAR_REAL2_CDF to read "real"
!                        scalars, as needed for WRFv2.2.  Added new routine
!                        GET_DIM_INT_CDF to retrieve netCDF dimensions.
!                        Changed internal error handling so that errors are
!                        passed back using a non-zero RCODE for dimension
!                        mismatches in addition to netCDF errors.  (T. Otte)
!           12 Feb 2010  Removed unused variable ID_TIME from subroutine
!                        GET_GL_ATT_INT_CDF.  (T. Otte)
!           19 Mar 2010  Removed routines GET_DIMS_CDF, GET_DIM_INT_CDF,
!                        GET_GL_ATT_INT_CDF, GET_GL_ATT_REAL_CDF,
!                        GET_GL_ATT_TEXT_CDF, and GET_DIM_ATT_INT_CDF.
!                        Removed file open and close functions from all
!                        remaining routines, and changed input argument
!                        from FILENAME to CDFID.  (T. Otte)
!           31 Aug 2011  Updated netCDF to F90.  Removed GET_TIMES_CDF.
!                        Changed F77 character declarations to F90 standard.
!                        (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           02 Feb 2018  Added new routine GET_VAR_3D_INT_CDF.  (T. Spero)
!           22 Jun 2018  Changed module name from WRF_NETCDF to NETCDF_IO.
!                        (T. Spero)
!           19 Dec 2022  Brought module to the Canopy-App (P.C. Campbell)
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!    SUBROUTINE get_var_3d_real_cdf (cdfid, var, dum3d, it, rcode)
!
!        USE canopy_const_mod, ONLY: rk
!        USE netcdf
!
!        IMPLICIT NONE
!
!        INTEGER,           INTENT(IN)    :: cdfid
!        REAL(rk),          INTENT(OUT)   :: dum3d    ( : , : , : )
!        INTEGER                          :: id_data
!        INTEGER,           INTENT(IN)    :: it
!        INTEGER                          :: nx
!        INTEGER                          :: ny
!        INTEGER                          :: nz
!        INTEGER,           INTENT(OUT)   :: rcode
!        CHARACTER(LEN=*),  INTENT(IN)    :: var
!
!
!        nx = SIZE(dum3d,2)
!        ny = SIZE(dum3d,1)
!        nz = SIZE(dum3d,3)
!
!        rcode = nf90_inq_varid (cdfid, var, id_data)
!        IF ( rcode /= nf90_noerr ) RETURN
!
!        rcode = nf90_get_var (cdfid, id_data, dum3d, start=(/1,1,1,it/),  &
!            count=(/nx,ny,nz,1/))
!
!        IF ( rcode /= nf90_noerr ) then
!            print*,'read error ',cdfid,var
!            print*,'nx,ny,nz=',nx,ny,nz
!        endif
!    END SUBROUTINE get_var_3d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!    SUBROUTINE get_var_3d_int_cdf (cdfid, var, idum3d, it, rcode)
!
!        USE netcdf
!
!        IMPLICIT NONE
!
!        INTEGER,           INTENT(IN)    :: cdfid
!        INTEGER                          :: id_data
!        INTEGER,           INTENT(OUT)   :: idum3d    ( : , : , : )
!        INTEGER,           INTENT(IN)    :: it
!        INTEGER                          :: nx
!        INTEGER                          :: ny
!        INTEGER                          :: nz
!        INTEGER,           INTENT(OUT)   :: rcode
!        CHARACTER(LEN=*),  INTENT(IN)    :: var
!
!        nx = SIZE(idum3d,2)
!        ny = SIZE(idum3d,1)
!        nz = SIZE(idum3d,3)
!
!        rcode = nf90_inq_varid (cdfid, var, id_data)
!        IF ( rcode /= nf90_noerr ) RETURN
!
!        rcode = nf90_get_var (cdfid, id_data, idum3d, start=(/1,1,1,it/),  &
!            count=(/nx,ny,nz,1/))
!
!    END SUBROUTINE get_var_3d_int_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_2d_real_cdf (cdfid, var, dum2d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum2d    ( : , : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(dum2d,1)
        ny = SIZE(dum2d,2)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) then
            print*,'can not find variable ',trim(var)
            RETURN
        endif
        rcode = nf90_get_var (cdfid, id_data, dum2d, start=(/1,1,it/),  &
            count=(/nx,ny,1/))
        IF ( rcode /= nf90_noerr ) print*,'read error for ',trim(var),&
            ' id_data,it,nx,ny=',id_data,it,nx,ny

    END SUBROUTINE get_var_2d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_2d_int_cdf (cdfid, var, idum2d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: idum2d   ( : , : )
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(idum2d,1)
        ny = SIZE(idum2d,2)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, idum2d, start=(/1,1,it/),  &
            count=(/nx,ny,1/))

    END SUBROUTINE get_var_2d_int_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_1d_real_cdf (cdfid, var, dum1d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum1d    ( : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(dum1d)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, dum1d, start=(/1,it/),  &
            count=(/nx,1/))

    END SUBROUTINE get_var_1d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_1d_int_cdf (cdfid, var, idum1d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: idum1d   ( : )
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(idum1d)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, idum1d, start=(/1,it/),  &
            count=(/nx,1/))

    END SUBROUTINE get_var_1d_int_cdf

!-------
!    SUBROUTINE get_var_1d_double_cdf (cdfid, var, dum1d, it, rcode)
!
!        USE netcdf
!
!        IMPLICIT NONE
!
!        INTEGER,           INTENT(IN)    :: cdfid
!        REAL,  INTENT(OUT)               :: dum1d    ( : )
!        INTEGER                          :: id_data
!        INTEGER,           INTENT(IN)    :: it
!        INTEGER                          :: nx
!        INTEGER,           INTENT(OUT)   :: rcode
!        CHARACTER(LEN=*),  INTENT(IN)    :: var
!        double precision, allocatable  :: dbtmp(:)
!
!        nx = SIZE(dum1d)
!        allocate(dbtmp(nx))
!
!        rcode = nf90_inq_varid (cdfid, var, id_data)
!        IF ( rcode /= nf90_noerr ) RETURN
!
!        rcode = nf90_get_var (cdfid, id_data, dbtmp, start=(/1,it/),  &
!            count=(/nx,1/))
!        dum1d(:)=sngl(dbtmp(:))
!    END SUBROUTINE get_var_1d_double_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!    SUBROUTINE get_var_real_cdf (cdfid, var, scalar, rcode)
!
!        USE canopy_const_mod, ONLY: rk
!        USE netcdf
!
!        IMPLICIT NONE
!
!        INTEGER,           INTENT(IN)    :: cdfid
!        INTEGER                          :: id_data
!        INTEGER,           INTENT(OUT)   :: rcode
!        REAL(rk),          INTENT(OUT)   :: scalar
!        CHARACTER(LEN=*),  INTENT(IN)    :: var
!
!        rcode = nf90_inq_varid (cdfid, var, id_data)
!        IF ( rcode /= nf90_noerr ) RETURN
!
!        rcode = nf90_get_var (cdfid, id_data, scalar)
!
!    END SUBROUTINE get_var_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


    SUBROUTINE canopy_outncf_init

        !-------------------------------------------------------------------------------
        ! Name:     Initialize output NetCDF data structures.
        ! Purpose:  Initializes output NetCDF structures.
        ! Revised:  21 Dec 2022  Initial version.  (P. C. Campbell)
        !-------------------------------------------------------------------------------

        use canopy_const_mod, ONLY: fillreal      !constants for canopy models
        use canopy_coord_mod   !main canopy coordinate descriptions
        use canopy_canvars_mod !main canopy variables descriptions

        IMPLICIT NONE

        !-------------------------------------------------------------------------------
        ! Time-independent 2d fields at cell centers.
        !-------------------------------------------------------------------------------

        g_lat%fld = fillreal
        g_lat%fldname = 'lat'
        g_lat%long_name = 'latitude at cell centers'
        g_lat%units = 'degree_north'
        g_lat%dimnames(1) = 'nlon'
        g_lat%dimnames(2) = 'nlat'
        g_lat%istart(1) = 1
        g_lat%istart(2) = 1
        g_lat%iend(1) = nlon
        g_lat%iend(2) = nlat

        g_lon%fld = fillreal
        g_lon%fldname = 'lon'
        g_lon%long_name = 'longitude at cell centers'
        g_lon%units = 'degree_east'
        g_lon%dimnames(1) = 'nlon'
        g_lon%dimnames(2) = 'nlat'
        g_lon%istart(1) = 1
        g_lon%istart(2) = 1
        g_lon%iend(1) = nlon
        g_lon%iend(2) = nlat

        !-------------------------------------------------------------------------------
        ! Time-varying 2d fields at cell centers.
        !-------------------------------------------------------------------------------

        c_waf%fld = fillreal
        c_waf%fldname = 'waf'
        c_waf%long_name = 'wind adjustment factor'
        c_waf%units = '1'
        c_waf%fillvalue = fillreal
        c_waf%dimnames(1) = 'nlon'
        c_waf%dimnames(2) = 'nlat'
        c_waf%istart(1) = 1
        c_waf%istart(2) = 1
        c_waf%iend(1) = nlon
        c_waf%iend(2) = nlat

        c_flameh%fld = fillreal
        c_flameh%fldname = 'flameh'
        c_flameh%long_name = 'flame height'
        c_flameh%units = 'm'
        c_flameh%fillvalue = fillreal
        c_flameh%dimnames(1) = 'nlon'
        c_flameh%dimnames(2) = 'nlat'
        c_flameh%istart(1) = 1
        c_flameh%istart(2) = 1
        c_flameh%iend(1) = nlon
        c_flameh%iend(2) = nlat


        !-------------------------------------------------------------------------------
        ! Time-varying 3d fields at cell centers.
        !-------------------------------------------------------------------------------

        c_canwind%fld = fillreal
        c_canwind%fldname = 'canwind'
        c_canwind%long_name = 'above/below canopy wind speed'
        c_canwind%units = 'm s-1'
        c_canwind%fillvalue = fillreal
        c_canwind%dimnames(1) = 'nlon'
        c_canwind%dimnames(2) = 'nlat'
        c_canwind%dimnames(3) = 'nlays'
        c_canwind%istart(1) = 1
        c_canwind%istart(2) = 1
        c_canwind%istart(3) = 1
        c_canwind%iend(1) = nlon
        c_canwind%iend(2) = nlat
        c_canwind%iend(3) = modlays

        c_Kz%fld = fillreal
        c_Kz%fldname = 'kz'
        c_Kz%long_name = 'eddy diffusivity'
        c_Kz%units = 'm2 s-1'
        c_Kz%fillvalue = fillreal
        c_Kz%dimnames(1) = 'nlon'
        c_Kz%dimnames(2) = 'nlat'
        c_Kz%dimnames(3) = 'nlays'
        c_Kz%istart(1) = 1
        c_Kz%istart(2) = 1
        c_Kz%istart(3) = 1
        c_Kz%iend(1) = nlon
        c_Kz%iend(2) = nlat
        c_Kz%iend(3) = modlays

        c_rjcf%fld = fillreal
        c_rjcf%fldname = 'rjcf'
        c_rjcf%long_name = 'photolysis attenuation correction factor'
        c_rjcf%units = '1'
        c_rjcf%fillvalue = fillreal
        c_rjcf%dimnames(1) = 'nlon'
        c_rjcf%dimnames(2) = 'nlat'
        c_rjcf%dimnames(3) = 'nlays'
        c_rjcf%istart(1) = 1
        c_rjcf%istart(2) = 1
        c_rjcf%istart(3) = 1
        c_rjcf%iend(1) = nlon
        c_rjcf%iend(2) = nlat
        c_rjcf%iend(3) = modlays

    END SUBROUTINE canopy_outncf_init


    SUBROUTINE canopy_outncf_alloc

        !-------------------------------------------------------------------------------
        ! Name:     Allocate Arrays for NetCDF output Dimensions
        ! Purpose:  Allocate arrays for NetCDF output Dimensions
        ! Revised:  23 Dec 2022  Original version.  (P. C. Campbell)
        !-------------------------------------------------------------------------------

        USE canopy_coord_mod
        USE canopy_canvars_mod

        IMPLICIT NONE

        INTEGER                      :: nn

        !-------------------------------------------------------------------------------
        ! Time-independent 2d fields at cell centers.
        !-------------------------------------------------------------------------------

        nfld2dxy = 0

        nfld2dxy = nfld2dxy + 2   !LAT, LON

        ALLOCATE ( fld2dxy ( nfld2dxy ) )

        DO nn = 1, nfld2dxy
            ALLOCATE ( fld2dxy(nn)%fld(nlon,nlat) )
        ENDDO

        g_lat     => fld2dxy( 1)
        g_lon     => fld2dxy( 2)

        !-------------------------------------------------------------------------------
        ! Time-varying 2d fields at cell centers.
        !-------------------------------------------------------------------------------

        nfld2dxyt = 0

        nfld2dxyt = nfld2dxyt + 1  !WAF

        nfld2dxyt = nfld2dxyt + 1  !FLAMEH

        ALLOCATE ( fld2dxyt ( nfld2dxyt ) )

        DO nn = 1, nfld2dxyt
            ALLOCATE ( fld2dxyt(nn)%fld(nlon,nlat) )
        ENDDO

        c_waf       => fld2dxyt( 1 )
        c_flameh    => fld2dxyt( 2 )

        !-------------------------------------------------------------------------------
        ! Time-varying 3d fields at cell centers.
        !-------------------------------------------------------------------------------

        nfld3dxyzt = 0

        nfld3dxyzt = nfld3dxyzt + 1 !canwind

        nfld3dxyzt = nfld3dxyzt + 1 !KZ

        nfld3dxyzt = nfld3dxyzt + 1 !RJCF

        ALLOCATE ( fld3dxyzt ( nfld3dxyzt ) )

        DO nn = 1, nfld3dxyzt
            ALLOCATE ( fld3dxyzt(nn)%fld(nlon,nlat,modlays) )
        ENDDO

        c_canwind    => fld3dxyzt( 1 )
        c_Kz         => fld3dxyzt( 2 )
        c_rjcf       => fld3dxyzt( 3 )

    END SUBROUTINE canopy_outncf_alloc

    SUBROUTINE canopy_outncfglobal (cdfid_in, fl)

        !-------------------------------------------------------------------------------
        ! Name:     Output netCDF Global Attributes
        ! Purpose:  Write netCDF global attributes.
        ! Revised:  23 Dec 2022  Original version.  (P. C. Campbell)
        !-------------------------------------------------------------------------------

        USE canopy_coord_mod
        USE canopy_canvars_mod
        USE netcdf

        IMPLICIT NONE

        INTEGER,            INTENT(IN)  :: cdfid_in
        !  CHARACTER(LEN=32)               :: cstr
        CHARACTER(LEN=256), INTENT(IN)  :: fl
        CHARACTER(LEN=32),  PARAMETER   :: pname      = 'CANOPY_OUTNCFGLOBAL'
        INTEGER                         :: rcode
        CHARACTER(LEN=32)               :: var

        !-------------------------------------------------------------------------------
        ! Error, warning, and informational messages.
        !-------------------------------------------------------------------------------

        CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR CREATING ATTRIBUTE FOR', a, &
        & /, 1x, '***   IN FILE ', a, &
        & /, 1x, '***   ', a, &
        & /, 1x, 70('*'))"

        !-------------------------------------------------------------------------------
        ! Define global attributes.
        !-------------------------------------------------------------------------------

        var = "PROGNAME"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, progname)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "VERSION"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, ver)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "CODE_DATE"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, vdate)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        !  var = "INPUT_MODEL"
        !  IF ( met_model == 2 ) THEN
        !    cstr = "WRF ARW " // TRIM(met_release)
        !  ELSE IF ( met_model == 3 ) THEN
        !    cstr = "FV3 " // TRIM(met_release)
        !  ELSE
        !    cstr = " "
        !  ENDIF
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, cstr)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF

        var = "im"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, nlon)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "jm"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, nlat)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        var = "levels"
        rcode = nf90_put_att (cdfid_in, nf90_global, var, modlays)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
                TRIM(nf90_strerror(rcode))
            CALL exit (2)
        ENDIF

        !  var = "NTHIK"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, nthik)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "GDTYP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, gdtyp_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "P_ALP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, p_alp_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "P_BET"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, p_bet_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "P_GAM"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, p_gam_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "XCENT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, xcent_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "YCENT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, ycent_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "XORIG"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, xorig_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "YORIG"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, yorig_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "XCELL"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, xcell_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "YCELL"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, ycell_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGTYP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vgtyp_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGTOP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vgtop_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGTOP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vgtop_gd)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "VGLVLS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, vglvs_gd(:))
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "CEN_LAT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_cen_lat)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "CEN_LON"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_cen_lon)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "TRUELAT1"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_tru1)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "TRUELAT2"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_tru2)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MOAD_CEN_LAT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_proj_clat)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_REF_LAT"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_ref_lat)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "STAND_LON"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_proj_clon)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "DX"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_resoln)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "DY"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_resoln)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "PTOP"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_ptop)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_CUMULUS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_cumulus)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_SHAL_CU"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_shal_cu)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_MICROPHYS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_expl_moist)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_LW_RAD"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_lw_rad)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_SW_RAD"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_sw_rad)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_PBL"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_pbl)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_SFC_LAY"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_sfc_lay)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_LSM"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_soil_lsm)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_URBAN_PHYS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_urban_phys)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "LAND_USE_SOURCE"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_lu_src)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "MET_FDDA_3DAN"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_3dan)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  IF ( met_fdda_3dan > 0 ) THEN  ! 3d nudging (any variety)
        !
        !    var = "MET_FDDA_3DAN_WIND"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gv3d)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_3DAN_TEMP"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gt3d)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_3DAN_MOIS"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gq3d)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    IF ( met_fdda_3dan == 2 ) THEN  ! spectral nudging only
        !      var = "MET_FDDA_3DAN_GEOP"
        !      rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gph3d)
        !      IF ( rcode /= nf90_noerr ) THEN
        !        WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                        TRIM(nf90_strerror(rcode))
        !        CALL graceful_stop (pname)
        !      ENDIF
        !    ENDIF
        !
        !  ENDIF  ! 3d nudging
        !
        !  var = "MET_FDDA_SFAN"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_sfan)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  IF ( met_fdda_sfan > 0 ) THEN  ! surface nudging
        !
        !    var = "MET_FDDA_SFAN_WIND"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gvsfc)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_SFAN_TEMP"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gtsfc)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_SFAN_MOIS"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_gqsfc)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !  ENDIF  ! surface nudging
        !
        !  var = "MET_FDDA_OBS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_obs)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  IF ( met_fdda_obs > 0 ) THEN  ! observation nudging
        !
        !    var = "MET_FDDA_OBS_WIND"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_giv)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_OBS_TEMP"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_git)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !    var = "MET_FDDA_OBS_MOIS"
        !    rcode = nf90_put_att (cdfid_in, nf90_global, var, met_fdda_giq)
        !    IF ( rcode /= nf90_noerr ) THEN
        !      WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                      TRIM(nf90_strerror(rcode))
        !      CALL graceful_stop (pname)
        !    ENDIF
        !
        !  ENDIF  ! surface nudging
        !
        !  var = "MET_HYBRID"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, met_hybrid)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF
        !
        !  var = "EARTH_RADIUS"
        !  rcode = nf90_put_att (cdfid_in, nf90_global, var, eradm)
        !  IF ( rcode /= nf90_noerr ) THEN
        !    WRITE (6,f9300) TRIM(pname), TRIM(var), TRIM(fl),  &
        !                    TRIM(nf90_strerror(rcode))
        !    CALL graceful_stop (pname)
        !  ENDIF

    END SUBROUTINE canopy_outncfglobal


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

            var = "grid_xt"
            rcode = nf90_def_dim (cdfid_m, TRIM(var), nlon, dim_nx)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (6,f9100) TRIM(pname), TRIM(var), TRIM(fl),  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            var = "grid_yt"
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
                rcode = nf90_put_att (cdfid_m, id_fld(nn), '_FillValue',  &
                    fld2dxyt(n)%fillvalue)
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
                rcode = nf90_put_att (cdfid_m, id_fld(nn), '_FillValue',  &
                    fld3dxyzt(n)%fillvalue)
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
            if (ifcanwind .or. ifcanwaf) then
                c_waf%fld = waf_2d
                c_flameh%fld = flameh_2d
            end if

            !-------------------------------------------------------------------------------
            ! Time-varying 3d fields at cell centers.
            !-------------------------------------------------------------------------------
            if (ifcanwind .or. ifcanwaf) then
                c_canwind%fld = canWIND_3d
            end if
            if (ifcaneddy) then
                c_Kz%fld      = Kz_3d
            end if
            if (ifcanphot) then
                c_rjcf%fld    = rjcf_3d
            end if

            !-------------------------------------------------------------------------------
            ! Write variables.
            !-------------------------------------------------------------------------------

            write(*,*)  'Writing NetCDF Output'
            write(*,*)  '-------------------------------'

            !-------------------------------------------------------------------------------
            ! Time-independent 2d fields at cell centers.
            !-------------------------------------------------------------------------------
            write(*,*)  'Writing Time-independent 2d fields'
            write(*,*)  '-------------------------------'

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
            write(*,*)  'Writing Time-varying 2d fields'
            write(*,*)  '-------------------------------'

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
            ! Time-varying 3d fields at cell centers.
            !-------------------------------------------------------------------------------
            write(*,*)  'Writing Time-varying 3d fields'
            write(*,*)  '-------------------------------'

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


    SUBROUTINE canopy_read_ncf(infile)

        !-------------------------------------------------------------------------------
        ! Name:     Read Canopy Met/Sfc Inputs from NetCDF
        ! Purpose:  Read Canopy Met/Sfc Inputs from NetCDF
        ! Revised:  19 Dec 2022  Original version.  (P.C. Campbell)
        !-------------------------------------------------------------------------------
        USE canopy_canopts_mod !main canopy option descriptions
        USE canopy_coord_mod
        USE canopy_canmet_mod
        USE netcdf

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT( IN )   :: infile
        CHARACTER(LEN=*), PARAMETER      :: pname      = 'canopy_read_ncf'
        CHARACTER(LEN=24)                :: rddate
        CHARACTER(LEN=32)                :: date_init
        INTEGER                          :: cdfid,rcode,varid,it
        double precision                 :: rdtime

        !-------------------------------------------------------------------------------
        ! Error, warning, and informational messages.
        !-------------------------------------------------------------------------------

        !    CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
        !    & /, 1x, '*** SUBROUTINE: ', a, &
        !    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM NCF INPUT FILE', &
        !    & /, 1x, '***   VARIABLE = ', a, &
        !    & /, 1x, '***   RCODE = ', a, &
        !    & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9410 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR RETRIEVING NCF ID FROM NCF INPUT FILE', &
        & /, 1x, '***   VARIABLE = ', a, &
        & /, 1x, '***   NCF: ', a, &
        & /, 1x, 70('*'))"

        !    CHARACTER(LEN=256), PARAMETER :: f9420 = "(/, 1x, 70('*'), &
        !    & /, 1x, '*** SUBROUTINE: ', a, &
        !    & /, 1x, '***   ERROR INQUIRING ABOUT VAR IN NCF INPUT FILE', &
        !    & /, 1x, '***   VARIABLE = ', a, &
        !    & /, 1x, '***   NCF: ', a, &
        !    & /, 1x, 70('*'))"

        !    CHARACTER(LEN=256), PARAMETER :: f9430 = "(/, 1x, 70('*'), &
        !    & /, 1x, '*** SUBROUTINE: ', a, &
        !    & /, 1x, '***   ERROR RETRIEVING DIMS FROM NCF INPUT FILE', &
        !    & /, 1x, '***   VARIABLE = ', a, &
        !    & /, 1x, '***   NCF: ', a, &
        !    & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9440 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR OPENING NCF INPUT FILE', &
        & /, 1x, '***   FILE = ', a, &
        & /, 1x, '***   NCF: ', a, &
        & /, 1x, 70('*'))"

        CHARACTER(LEN=256), PARAMETER :: f9450 = "(/, 1x, 70('*'), &
        & /, 1x, '*** SUBROUTINE: ', a, &
        & /, 1x, '***   ERROR RETRIEVING VAR ATT IN NCF INPUT FILE', &
        & /, 1x, '***   VARIABLE = ', a, &
        & /, 1x, '***   NCF: ', a, &
        & /, 1x, 70('*'))"

        !-------------------------------------------------------------------------------
        ! Open NCF input file
        !-------------------------------------------------------------------------------

        rcode = nf90_open (infile, nf90_nowrite, cdfid)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9440) TRIM(pname), TRIM(infile)
            call exit(2)
        ENDIF

        !-------------------------------------------------------------------------------
        ! Extract and check input lat/lon dimensions and time information.
        !-------------------------------------------------------------------------------
        !time information
        rcode = nf90_inq_varid (cdfid, 'time', varid)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9410) TRIM(pname), 'time',  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        rcode = nf90_get_att (cdfid, varid, 'units', date_init)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9450) TRIM(pname), 'SIMULATION_START_DATE',  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        rddate = date_init(13:31) // '.0000'
        rddate(11:11)='-'

        rcode=nf90_get_var(cdfid,varid,rdtime)
        IF ( rcode /= nf90_noerr ) THEN
            WRITE (*,f9410) TRIM(pname), 'time',  &
                TRIM(nf90_strerror(rcode))
            CALL exit(2)
        ENDIF

        !grid 1D or 2D  lat/lon variables
        it = int(rdtime) + 1

        if (infmt_opt .eq. 0) then !Input format is 2D

            CALL get_var_2d_real_cdf (cdfid, 'lat', variables_2d%lat, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lon',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%lat=reshape(variables_2d%lat,[size(variables_2d%lat)])
            CALL get_var_2d_real_cdf (cdfid, 'lon', variables_2d%lon, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lon',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%lon=reshape(variables_2d%lon,[size(variables_2d%lon)])
            !Canopy input met/sfc variables
            !Clumping index
            CALL get_var_2d_real_cdf (cdfid, 'clu', variables_2d%clu, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'clu',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%clu=reshape(variables_2d%clu,[size(variables_2d%clu)])
            !Cosine of solar zenith angle
            CALL get_var_2d_real_cdf (cdfid, 'csz', variables_2d%csz, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'csz',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%csz=reshape(variables_2d%csz,[size(variables_2d%csz)])
            !Forest Fraction
            CALL get_var_2d_real_cdf (cdfid, 'ffrac', variables_2d%ffrac, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ffrac',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%ffrac=reshape(variables_2d%ffrac,[size(variables_2d%ffrac)])
            !Forest canopy height
            CALL get_var_2d_real_cdf (cdfid, 'fh', variables_2d%fh, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'fh',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%fh=reshape(variables_2d%fh,[size(variables_2d%fh)])
            !Fire Radiative Power
            CALL get_var_2d_real_cdf (cdfid, 'frp', variables_2d%frp, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'frp',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%frp=reshape(variables_2d%frp,[size(variables_2d%frp)])
            !Reference height above canopy
            CALL get_var_2d_real_cdf (cdfid, 'href', variables_2d%href, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'href',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%href=reshape(variables_2d%href,[size(variables_2d%href)])
            !Leaf Area Index
            CALL get_var_2d_real_cdf (cdfid, 'lai', variables_2d%lai, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lai',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%lai=reshape(variables_2d%lai,[size(variables_2d%lai)])
            !Monin-Obukhov Length
            CALL get_var_2d_real_cdf (cdfid, 'mol', variables_2d%mol, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'mol',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%mol=reshape(variables_2d%mol,[size(variables_2d%mol)])
            !Friction velocity
            CALL get_var_2d_real_cdf (cdfid, 'fricv', variables_2d%fricv, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'fricv',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%fricv=reshape(variables_2d%fricv,[size(variables_2d%fricv)])
            !Reference U Wind Speed (at HREF)
            CALL get_var_2d_real_cdf (cdfid, 'ugrd10m', variables_2d%ugrd10m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ugrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%ugrd10m=reshape(variables_2d%ugrd10m,[size(variables_2d%ugrd10m)])
            !Reference V Wind Speed (at HREF)
            CALL get_var_2d_real_cdf (cdfid, 'vgrd10m', variables_2d%vgrd10m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vgrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%vgrd10m=reshape(variables_2d%vgrd10m,[size(variables_2d%vgrd10m)])
            !Surface (veg+soil) Roughness Length
            CALL get_var_2d_real_cdf (cdfid, 'sfcr', variables_2d%sfcr, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'sfcr',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%sfcr=reshape(variables_2d%sfcr,[size(variables_2d%sfcr)])
            !Vegetation Type
            CALL get_var_2d_int_cdf (cdfid, 'vtype', variables_2d%vtype, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vtype',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Also reshape to 1D array for 1D calculation and output
            variables%vtype=reshape(variables_2d%vtype,[size(variables_2d%vtype)])

        else if (infmt_opt .eq. 1) then !Input format is 1D

            CALL get_var_1d_real_cdf (cdfid, 'lat', variables%lat, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lat',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            CALL get_var_1d_real_cdf (cdfid, 'lon', variables%lon, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lon',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

            !Canopy input met/sfc variables
            !Clumping index
            CALL get_var_1d_real_cdf (cdfid, 'clu', variables%clu, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'clu',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Cosine of solar zenith angle
            CALL get_var_1d_real_cdf (cdfid, 'csz', variables%csz, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'csz',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Forest Fraction
            CALL get_var_1d_real_cdf (cdfid, 'ffrac', variables%ffrac, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ffrac',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Forest canopy height
            CALL get_var_1d_real_cdf (cdfid, 'fh', variables%fh, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'fh',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Fire Radiative Power
            CALL get_var_1d_real_cdf (cdfid, 'frp', variables%frp, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'frp',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Reference height above canopy
            CALL get_var_1d_real_cdf (cdfid, 'href', variables%href, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'href',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Leaf Area Index
            CALL get_var_1d_real_cdf (cdfid, 'lai', variables%lai, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'lai',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Monin-Obukhov Length
            CALL get_var_1d_real_cdf (cdfid, 'mol', variables%mol, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'mol',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Friction velocity
            CALL get_var_1d_real_cdf (cdfid, 'fricv', variables%fricv, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'fricv',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Reference U Wind Speed (at HREF)
            CALL get_var_1d_real_cdf (cdfid, 'ugrd10m', variables%ugrd10m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'ugrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Reference V Wind Speed (at HREF)
            CALL get_var_1d_real_cdf (cdfid, 'vgrd10m', variables%vgrd10m, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vgrd10m',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Surface (veg+soil) Roughness Length
            CALL get_var_1d_real_cdf (cdfid, 'sfcr', variables%sfcr, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'sfcr',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF
            !Vegetation Type
            CALL get_var_1d_int_cdf (cdfid, 'vtype', variables%vtype, it, rcode)
            IF ( rcode /= nf90_noerr ) THEN
                WRITE (*,f9410) TRIM(pname), 'vtype',  &
                    TRIM(nf90_strerror(rcode))
                CALL exit(2)
            ENDIF

        else
            write(*,*)  'Wrong INFMT_OPT choice of ', infmt_opt, ' in namelist...exiting'
            call exit(2)
        end if !Input Format (1D or 2D)

    END SUBROUTINE canopy_read_ncf

END MODULE canopy_ncf_io_mod
