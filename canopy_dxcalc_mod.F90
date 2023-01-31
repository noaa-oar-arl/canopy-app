module canopy_dxcalc_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_CALCDX(DXOPT, DXSET, NLAT, NLON, LAT, LON, DX )

!-----------------------------------------------------------------------

! Description:
!     computes great circle distance or the orthodromic distance using Haversine
!     formula

! Preconditions:
!     user dx_opt, dx_set, nlat, nlon, lat, and lon

! Subroutines and Functions Called:

! Revision History:
!     Prototype 10/22 by PCC
!     Oct 2022 P.C. Campbell: Initial version
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk       !constants for canopy models
        use canopy_utils_mod, ONLY: CalcDX   !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,    INTENT( IN )  :: DXOPT           ! User DX calculation option
        REAL(RK),   INTENT( IN )  :: DXSET           ! User DX set value if cannot calculate (m)
        INTEGER,    INTENT( IN )  :: NLAT            ! Number of latitude grid cells/points
        INTEGER,    INTENT( IN )  :: NLON            ! Number of longitude grid cells/points
        REAL(RK),   INTENT( IN )  :: LAT(:)          ! Model Latitudes (degrees)
        REAL(RK),   INTENT( IN )  :: LON(:)          ! Model Longitudes (degrees)
        REAL(RK),   INTENT( OUT ) :: DX(:)           ! Distance between two points (m)

!     Local variables
        integer  ::    loc                           ! NLAT*NLON

        do loc=1, NLAT*NLON

            if (DXOPT .eq. 0) then !user set to calculate dx grid cell distance from grid lons
                if (NLON .gt. 1) then !convert grid points to distances using Haversine formula (m)
                    if (loc .lt. NLAT*NLON) then !inside domain
                        DX(loc) = CalcDX(LAT(loc), abs(LON(loc+1) - LON(loc)))
                    else !at the domain edge --set to loc-1
                        DX(loc) = DX(NLAT*NLON-1)
                    end if
                else                  !single grid cell/point, use namelist defined dx resolution (m) for cell
                    write(*,*)  'DX_OPT set to calc, but nlon or nlat  <= 1...setting dx = ', &
                        DXSET, ' from namelist'
                    DX(loc) = DXSET
                end if
            else ! user set dx_set from namelist
                DX(loc) = DXSET
            end if
        end do

    END SUBROUTINE CANOPY_CALCDX
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_CALCDX_2D(DXOPT, DXSET, NLAT, NLON, LAT, LON, DX )

!-----------------------------------------------------------------------

! Description:
!     computes great circle distance or the orthodromic distance using Haversine
!     formula

! Preconditions:
!     user dx_opt, dx_set, nlat, nlon, lat, and lon

! Subroutines and Functions Called:

! Revision History:
!     Prototype 10/22 by PCC
!     Oct 2022 P.C. Campbell: Initial version
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk       !constants for canopy models
        use canopy_utils_mod, ONLY: CalcDX   !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,    INTENT( IN )  :: DXOPT           ! User DX calculation option
        REAL(RK),   INTENT( IN )  :: DXSET           ! User DX set value if cannot calculate (m)
        INTEGER,    INTENT( IN )  :: NLAT            ! Number of latitude grid cells/points
        INTEGER,    INTENT( IN )  :: NLON            ! Number of longitude grid cells/points
        REAL(RK),   INTENT( IN )  :: LAT(:,:)          ! Model Latitudes (degrees)
        REAL(RK),   INTENT( IN )  :: LON(:,:)          ! Model Longitudes (degrees)
        REAL(RK),   INTENT( OUT ) :: DX(:,:)           ! Distance between two points (m)

!     Local variables
        integer  ::    i,j                           ! NLON,NLAT

        do i=1, NLON
            do j=1, NLAT

                if (DXOPT .eq. 0) then !user set to calculate dx grid cell distance from grid lons
                    if (NLON .gt. 1) then !convert grid points to distances using Haversine formula (m)
                        if (i .lt. NLON ) then !inside LON inside domain
                            DX(i,j) = CalcDX(LAT(i,j), abs(LON(i+1,j) - LON(i,j)))
                        else if (i .eq. NLON ) then !at the domain edge --set to NLON-1
                            DX(i,j) = DX(NLON-1,j)
                        else if (j .eq. NLAT ) then !at the domain edge --set to NLAT-1
                            DX(i,j) = DX(i,NLAT-1)
                        end if
                    else                  !single grid cell/point, use namelist defined dx resolution (m) for cell
                        write(*,*)  'DX_OPT_2D set to calc, but nlon or nlat  <= 1...setting dx = ', &
                            DXSET, ' from namelist'
                        DX(i,j) = DXSET
                    end if
                else ! user set dx_set from namelist
                    DX(i,j) = DXSET
                end if
            end do !LAT
        end do   !LON

    END SUBROUTINE CANOPY_CALCDX_2D
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module canopy_dxcalc_mod
