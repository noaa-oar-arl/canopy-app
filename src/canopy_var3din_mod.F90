module canopy_var3din_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_PAVD2FAFRAC ( MODRES, FCH, ZHC, &
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
        use canopy_const_mod, ONLY: rk                                           !constants for canopy models
        use canopy_utils_mod, ONLY: interp_linear1_internal,IntegrateTrapezoid   !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: MODRES                ! Canopy model input vertical resolution (m)
        REAL(RK),    INTENT( IN )  :: FCH                   ! Grid cell canopy height (m)  !From GEDI
        REAL(RK),    INTENT( IN )  :: ZHC(:)                ! z/h (dimensionless)
        REAL(RK),    INTENT( IN )  :: PAVD_IN(:)            ! Plant Area Volume Density (PAVD) profile (m2/m3)
        REAL(RK),    INTENT( IN )  :: PAVD_LEVS(:)          ! Associated mid-level heights associated with PAVD (m)
        REAL(RK),    INTENT( OUT ) :: FAFRACZINT(:)         ! integral of incremental fractional foliage shape function

!     Local variables
        INTEGER                     :: i, lev
        REAL(RK)                    :: PAI                  ! Plant Area Index, m2/m2
        REAL(RK)                    :: ZK(SIZE(ZHC))
        REAL(RK)                    :: PAVD_INTERP(SIZE(ZHC))

        REAL(RK)                    :: ZCANMAX               ! Height of maximum foliage area density (z/h) (nondimensional)
        REAL(RK)                    :: PAVD_BELOW            ! Average of PAVD below zcanmax
        REAL(RK)                    :: PAVD_ABOVE            ! Average of PAVD above zcanmax
        REAL(RK)                    :: SIGMAU                ! Standard deviation of shape function above zcanmax (z/h)
        REAL(RK)                    :: SIGMA1                ! Standard deviation of shape function below zcanmax (z/h)



!        REAL(RK), allocatable       :: fainc(:)              ! incremental foliage shape function
!        REAL(RK), allocatable       :: fafracz(:)            ! incremental fractional foliage shape function
!        REAL(RK)                    :: fatot                 ! integral of total fractional foliage shape function

        ZK = ZHC*FCH

!Interpolate the input PAVD at set levels to the user desired canopy model resolution
        PAVD_INTERP = 0.0_rk  !Initialize PAVD_INTERP = 0
        do lev=1, SIZE(PAVD_LEVS) - 1
            do i=2, SIZE(ZK)  !loop over only levels ABOVE ground
!                if (ZK(i) .le. FCH) then !constrain to less than the forest canopy height observed (e.g., GEDI CH)
                if (ZK(i) .le.  PAVD_LEVS(1)) then
                    PAVD_INTERP(i)   = PAVD_IN(1)
                end if
                if (ZK(i) .ge.  PAVD_LEVS(lev) .and. ZK(i) .le.  PAVD_LEVS(lev+1)) then
                    PAVD_INTERP(i)   = interp_linear1_internal((/ PAVD_LEVS(lev),PAVD_LEVS(lev+1) /), &
                        (/ PAVD_IN(lev),PAVD_IN(lev+1) /),ZK(i))
                end if
!                end if
            end do
        end do

!bottom-up loop allowing ZCANMAX, SIGMAU, and SIGMA1 calculation from observed PAVD
        do i=2, SIZE(ZK)
            if (PAVD_INTERP(i) .ge. maxval(PAVD_INTERP) ) then
                ZCANMAX = ZK(i)/FCH
                PAVD_BELOW = SUM(PAVD_INTERP(2:i)) / i
                SIGMA1 = (SQRT (SUM((PAVD_INTERP(2:i)-PAVD_BELOW)**2) / i))/FCH
                PAVD_ABOVE = SUM(PAVD_INTERP(i:SIZE(ZK))) / (SIZE(ZK)-i)
                SIGMAU = (SQRT (SUM((PAVD_INTERP(i:SIZE(ZK))-PAVD_ABOVE)**2) / (SIZE(ZK)-i)))/FCH
                exit
            end if
        end do
        print*, ZCANMAX,PAVD_BELOW,SIGMA1,PAVD_ABOVE,SIGMAU

!Integrate the PAVD_INTERP to get the total PAI
        PAI = IntegrateTrapezoid(ZK,PAVD_INTERP)

!Convert the PAVD_INTERP to FAFRACZINT
        FAFRACZINT = PAVD_INTERP * 0.0_rk !Initialize FAFRACZINT = 0
        do i=2, SIZE(ZK)  !loop over only levels ABOVE ground
            if (PAI .gt. 0.0_rk) then !check if PAI from GEDI PAVD > 0
                FAFRACZINT(i) = FAFRACZINT(i-1) + ((PAVD_INTERP(i)*MODRES)/PAI)
            end if
        end do

        !test debug prints...
!        print*, 'FCH = ', FCH
!        print*, 'pavd_orig = ', PAVD_IN
!        print*, 'pavd_levs = ', PAVD_LEVS
!        print*,'PAVD After Interpolation (canopy-app)-----'
!        print*,PAVD_INTERP
!        print*,'------------------------'
!        print*,'canopy-app levels-----'
!        print*,ZK
!        print*,'------------------------'
!        print*, 'PAI (Integrated PAVD) = ', PAI

    END SUBROUTINE CANOPY_PAVD2FAFRAC

end module canopy_var3din_mod
