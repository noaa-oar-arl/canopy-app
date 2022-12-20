
SUBROUTINE canopy_read_ncf(infile)

!-------------------------------------------------------------------------------
! Name:     Read Canopy Met/Sfc Inputs from NetCDF
! Purpose:  Read Canopy Met/Sfc Inputs from NetCDF
! Revised:  19 Dec 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------
    USE canopy_coord_mod
    USE canopy_canmet_mod
    USE canopy_ncf_io_mod  !main IO NetCDF reader
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

    !grid lat/lon variables
    it = int(rdtime) + 1

    CALL get_var_1d_real_cdf (cdfid, 'LAT', variables%lat, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'LAT',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF

    CALL get_var_1d_real_cdf (cdfid, 'LON', variables%lon, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'LON',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF

    !Canopy input met/sfc variables
    !Clumping index
    CALL get_var_1d_real_cdf (cdfid, 'CLU', variables%clu, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'CLU',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Cosine of solar zenith angle
    CALL get_var_1d_real_cdf (cdfid, 'CSZ', variables%csz, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'CSZ',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Forest Fraction
    CALL get_var_1d_real_cdf (cdfid, 'FFRAC', variables%ffrac, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'FFRAC',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Forest canopy height
    CALL get_var_1d_real_cdf (cdfid, 'FH', variables%fh, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'FH',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Fire Radiative Power
    CALL get_var_1d_real_cdf (cdfid, 'FRP', variables%frp, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'FRP',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Reference height above canopy
    CALL get_var_1d_real_cdf (cdfid, 'HREF', variables%href, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'HREF',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Leaf Area Index
    CALL get_var_1d_real_cdf (cdfid, 'LAI', variables%lai, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'LAI',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Monin-Obukhov Length
    CALL get_var_1d_real_cdf (cdfid, 'MOL', variables%mol, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'MOL',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Friction velocity
    CALL get_var_1d_real_cdf (cdfid, 'UST', variables%ust, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'UST',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Reference Wind Speed (at HREF)
    CALL get_var_1d_real_cdf (cdfid, 'WS', variables%ws, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'WS',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Surface (veg+soil) Roughness Length
    CALL get_var_1d_real_cdf (cdfid, 'Z0', variables%z0, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'Z0',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF
    !Vegetation Type
    CALL get_var_1d_int_cdf (cdfid, 'VTYPE', variables%vtype, it, rcode)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9410) TRIM(pname), 'VTYPE',  &
            TRIM(nf90_strerror(rcode))
        CALL exit(2)
    ENDIF

END SUBROUTINE canopy_read_ncf
