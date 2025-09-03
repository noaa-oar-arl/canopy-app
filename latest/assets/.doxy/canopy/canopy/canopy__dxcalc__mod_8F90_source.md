

# File canopy\_dxcalc\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_dxcalc\_mod.F90**](canopy__dxcalc__mod_8F90.md)

[Go to the documentation of this file](canopy__dxcalc__mod_8F90.md)


```Fortran



module canopy_dxcalc_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_calcdx(DXOPT, DXSET, NLAT, NLON, LAT, LON, DX )

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
        use canopy_utils_mod, ONLY: calcdx   !utilities for canopy models

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
        integer  ::    loc

        do loc=1, nlat*nlon

            if (dxopt .eq. 0) then !user set to calculate dx grid cell distance from grid lons
                if (nlon .gt. 1) then !convert grid points to distances using Haversine formula (m)
                    if (loc .lt. nlat*nlon) then !inside domain
                        dx(loc) = calcdx(lat(loc), abs(lon(loc+1) - lon(loc)))
                    else !at the domain edge --set to loc-1
                        dx(loc) = dx(nlat*nlon-1)
                    end if
                else                  !single grid cell/point, use namelist defined dx resolution (m) for cell
                    write(*,*)  'DX_OPT set to calc, but nlon or nlat  <= 1...setting dx = ', &
                        dxset, ' from namelist'
                    dx(loc) = dxset
                end if
            else ! user set dx_set from namelist
                dx(loc) = dxset
            end if
        end do

    END SUBROUTINE canopy_calcdx
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_calcdx_2d(DXOPT, DXSET, NLAT, NLON, LAT, LON, DX )

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
        use canopy_utils_mod, ONLY: calcdx   !utilities for canopy models

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
        integer  ::    i,j

        do i=1, nlon
            do j=1, nlat

                if (dxopt .eq. 0) then !user set to calculate dx grid cell distance from grid lons
                    if (nlon .gt. 1) then !convert grid points to distances using Haversine formula (m)
                        if (i .lt. nlon ) then !inside LON inside domain
                            dx(i,j) = calcdx(lat(i,j), abs(lon(i+1,j) - lon(i,j)))
                        else if (i .eq. nlon ) then !at the domain edge --set to NLON-1
                            dx(i,j) = dx(nlon-1,j)
                        else if (j .eq. nlat ) then !at the domain edge --set to NLAT-1
                            dx(i,j) = dx(i,nlat-1)
                        end if
                    else                  !single grid cell/point, use namelist defined dx resolution (m) for cell
                        write(*,*)  'DX_OPT_2D set to calc, but nlon or nlat  <= 1...setting dx = ', &
                            dxset, ' from namelist'
                        dx(i,j) = dxset
                    end if
                else ! user set dx_set from namelist
                    dx(i,j) = dxset
                end if
            end do !LAT
        end do   !LON

    END SUBROUTINE canopy_calcdx_2d
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module canopy_dxcalc_mod
```


