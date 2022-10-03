module canopy_calcdx_mod

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
                if (NLON .gt. 1 .and. NLON .gt. 1) then !convert grid points to distances using Haversine formula (m)
                    if (loc .lt. NLAT*NLON) then !inside domain
                        DX(loc) = CalcDX(LAT(loc),LAT(loc+1),LON(loc),LON(loc+1)) !Haversine Calc
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

end module canopy_calcdx_mod
