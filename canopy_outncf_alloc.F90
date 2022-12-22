
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

    ALLOCATE ( fld2dxyt ( nfld2dxyt ) )

    DO nn = 1, nfld2dxyt
        ALLOCATE ( fld2dxyt(nn)%fld(nlon,nlat) )
    ENDDO

    c_waf    => fld2dxyt( 1 )

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
