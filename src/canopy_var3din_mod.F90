module canopy_var3din_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_PAVD2FAFRAC ( MODLAYS, MODRES, FCH, ZHC, &
        PAVD_IN, PAVD_LEVS, FAFRACZINT )

!-----------------------------------------------------------------------

! Description:
!     computes integral of incremental fractional foliage shape function
!     from interpolated GEDI 3D PAVD profile.

! Preconditions:
!     height of canopy, PAVD, etc.

! Subroutines and Functions Called:

! Revision History:
!     Prototype 07/23 by PCC, based on GEDI PAVD
!     Jul 2023 P.C. Campbell: Initial version
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk                        !constants for canopy models
        use canopy_utils_mod, ONLY: interp_linear1_internal   !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )  :: MODLAYS               ! Number of total above and below canopy model layers
        REAL(RK),    INTENT( IN )  :: MODRES                ! Canopy model input vertical resolution (m)
        REAL(RK),    INTENT( IN )  :: FCH                   ! Grid cell canopy height (m)  !From GEDI
        REAL(RK),    INTENT( IN )  :: ZHC(:)                ! z/h (dimensionless)
        REAL(RK),    INTENT( IN )  :: PAVD_IN(:)            ! Plant Area Volume Density (PAVD) profile (m2/m3)
        REAL(RK),    INTENT( IN )  :: PAVD_LEVS(:)          ! Associated mid-level heights associated with PAVD (m)
        REAL(RK),    INTENT( OUT ) :: FAFRACZINT(:)         ! integral of incremental fractional foliage shape function

!     Local variables
        INTEGER                     :: i, lev
        REAL(RK)                    :: ZK(SIZE(ZHC))
        REAL(RK)                    :: PAVD_INTERP(SIZE(ZHC))
!        REAL(RK), allocatable       :: fainc(:)              ! incremental foliage shape function
!        REAL(RK), allocatable       :: fafracz(:)            ! incremental fractional foliage shape function
!        REAL(RK)                    :: fatot                 ! integral of total fractional foliage shape function

        ZK = ZHC*FCH

!Interpolate the input PAVD at set levels to the user desired canopy model resolution
        PAVD_INTERP = 0.0_rk  !Initialize PAVD_INTERP = 0
        do lev=1, SIZE(PAVD_LEVS) - 1
            do i=2, SIZE(ZK)  !loop over only levels ABOVE ground
                if (ZK(i) .le.  PAVD_LEVS(1)) then
                    PAVD_INTERP(i)   = PAVD_IN(1)
                end if
                if (ZK(i) .ge.  PAVD_LEVS(lev) .and. ZK(i) .le.  PAVD_LEVS(lev+1)) then
                    PAVD_INTERP(i)   = interp_linear1_internal((/ PAVD_LEVS(lev),PAVD_LEVS(lev+1) /), &
                        (/ PAVD_IN(lev),PAVD_IN(lev+1) /),ZK(i))
                end if
            end do
        end do

        FAFRACZINT = PAVD_INTERP  !TODO:  Convert PAVD_INTERP to FAFRACZINT

!        print*,'After interpolation-----'
!        print*,PAVD_INTERP
!        print*,'------------------------'

    END SUBROUTINE CANOPY_PAVD2FAFRAC

end module canopy_var3din_mod
