
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
    g_lat%fldname = 'LAT'
    g_lat%long_name = 'latitude at cell centers'
    g_lat%units = 'degrees_north'
    g_lat%dimnames(1) = 'nlon'
    g_lat%dimnames(2) = 'nlat'
    g_lat%istart(1) = 1
    g_lat%istart(2) = 1
    g_lat%iend(1) = nlon
    g_lat%iend(2) = nlat

    g_lon%fld = fillreal
    g_lon%fldname = 'LON'
    g_lon%long_name = 'longitude at cell centers'
    g_lon%units = 'degrees_east'
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
    c_waf%fldname = 'WAF'
    c_waf%long_name = 'wind adjustment factor'
    c_waf%units = '1'
    c_waf%fillvalue = fillreal
    c_waf%dimnames(1) = 'nlon'
    c_waf%dimnames(2) = 'nlat'
    c_waf%istart(1) = 1
    c_waf%istart(2) = 1
    c_waf%iend(1) = nlon
    c_waf%iend(2) = nlat

!-------------------------------------------------------------------------------
! Time-varying 3d fields at cell centers.
!-------------------------------------------------------------------------------

    c_canwind%fld = fillreal
    c_canwind%fldname = 'CANWIND'
    c_canwind%long_name = 'Above/below canopy wind speeds'
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
    c_Kz%fldname = 'KZ'
    c_Kz%long_name = 'eddy diffusivities'
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
    c_rjcf%fldname = 'RJCF'
    c_rjcf%long_name = 'photolysis attenuation correction factors'
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
