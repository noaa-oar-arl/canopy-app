module canopy_var3din_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_PAVD2FAFRAC ( ZCANMAX_IN, SIGMAU, SIGMA1, FCH, &
        ZHC, PAVD_IN, PAVD_LEVS, FAFRACZINT )


!-----------------------------------------------------------------------

! Description:
!     computes integral of incremental fractional foliage shape function
!     from interpolated GEDI 3D PAVD profile and ZCANMAX.

! Preconditions:
!     height of canopy, PAVD, etc.

! Subroutines and Functions Called:

! Revision History:
!     Prototype 07/23 by PCC, based on GEDI PAVD
!     Jul 2023 P.C. Campbell: Initial version
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk                                           !constants for canopy models
        use canopy_utils_mod, ONLY: interp_linear1_internal,IntegrateTrapezoid   !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: FCH                   ! Grid cell canopy height (m)  !From GEDI
        REAL(RK),    INTENT( IN )  :: ZCANMAX_IN            ! Input height of maximum foliage area density (z/h) (nondimensional)
        REAL(RK),    INTENT( IN )  :: SIGMAU                ! Standard deviation of shape function above zcanmax (z/h)
        REAL(RK),    INTENT( IN )  :: SIGMA1                ! Standard deviation of shape function below zcanmax (z/h)
        REAL(RK),    INTENT( IN )  :: ZHC(:)                ! z/h (dimensionless)
        REAL(RK),    INTENT( IN )  :: PAVD_IN(:)            ! Plant Area Volume Density (PAVD) profile (m2/m3)
        REAL(RK),    INTENT( IN )  :: PAVD_LEVS(:)          ! Associated mid-level heights associated with PAVD (m)
        REAL(RK),    INTENT( OUT ) :: FAFRACZINT(:)         ! integral of incremental fractional foliage shape function

!     Local variables
        INTEGER                     :: i, lev
        REAL(RK)                    :: ZK(SIZE(ZHC))
        REAL(RK)                    :: PAVD_INTERP(SIZE(ZHC))

        REAL(RK)                    :: ZCANMAX               ! Height of maximum foliage area density (z/h) (nondimensional)
        REAL(RK), allocatable       :: fainc(:)              ! incremental foliage shape function
        REAL(RK), allocatable       :: fafracz(:)            ! incremental fractional foliage shape function
        REAL(RK)                    :: fatot                 ! integral of total fractional foliage shape function

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

! ... initialize canopy profile dependent variables
        ZCANMAX             = 0.0_rk

!bottom-up loop allowing ZCANMAX calculation from observed PAVD
        do i=2, SIZE(ZK)
            if (PAVD_INTERP(i) .ge. maxval(PAVD_INTERP) ) then
                ZCANMAX = ZK(i)/FCH
                if(ZCANMAX .gt. 1.0_rk) then !if ZK (at max GEDI PAVD height) > GEDI FCH (inconsistent!)
                    ZCANMAX = ZCANMAX_IN       !override with Massman Input ZCANMAX
                end if
                exit
            end if
        end do

!Eqs. 1-2 on pg 595 of Massman et al., 2017

        if(.not.allocated(fainc))      allocate(fainc(SIZE(ZK)))
        if(.not.allocated(fafracz))    allocate(fafracz(SIZE(ZK)))

! ... initialize canopy profile dependent variables
        fainc             = 0.0_rk
        fafracz           = 0.0_rk

! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
        do i=1, SIZE(ZK)
            if (ZHC(i) >= ZCANMAX .and. ZHC(i) <= 1.0) then
                fainc(i) = exp((-1.0*((ZHC(i)-ZCANMAX)**2.0))/SIGMAU**2.0)
            else if (ZHC(i) >= 0.0 .and. ZHC(i) <= ZCANMAX) then
                fainc(i) = exp((-1.0*((ZCANMAX-ZHC(i))**2.0))/SIGMA1**2.0)
            end if
        end do
        fatot = IntegrateTrapezoid(ZHC,fainc)

! ... calculate plant distribution function
        do i=1, SIZE(ZK)
            fafracz(i) = fainc(i)/fatot
            FAFRACZINT(i) = IntegrateTrapezoid(ZHC(1:i),fafracz(1:i))
        end do

    END SUBROUTINE CANOPY_PAVD2FAFRAC

end module canopy_var3din_mod
