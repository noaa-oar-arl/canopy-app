module canopy_zpd_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_ZPD( ZHC, FCLAI, UBZREF, Z0GHC, &
        LAMDARS, CDRAG, PAI, d_h, zo_h )

!-----------------------------------------------------------------------

! Description:
!     computes zero plane displacement height and total surface roughness length basd on Massman et al. (2017).

! Preconditions:
!     foliage characteristics, drag, and in-canopy plant distribution functions

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone canopy wind model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk, vonk    !constants for canopy models
        use canopy_utils_mod     !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: ZHC     ( : )   ! SUB-CANOPY Total z/h layers (nondimensional)
        REAL(RK),    INTENT( IN )  :: FCLAI   ( : )   ! SUB-CANOPY Fractional (z) shapes of the
        ! plant surface distribution, i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )  :: UBZREF          ! Mean wind speed at zref-height of canopy top (m/s)
        REAL(RK),    INTENT( IN )  :: Z0GHC           ! Ratio of ground roughness length to canopy top height (nondimensional)
        REAL(RK),    INTENT( IN )  :: LAMDARS         ! Value representing influence of roughness sublayer (nondimensional)
        REAL(RK),    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
        REAL(RK),    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
        REAL(RK),    INTENT( OUT ) :: d_h             ! zero plane displacement height d/h (i.e., Product of dha*dhb)
        REAL(RK),    INTENT( OUT ) :: zo_h            ! Surface (soil+veg) roughness length, zo/h

!     Local variables
        real(rk)                   :: ustrmod         ! Friction Velocity parameterization (m/s)
        real(rk)                   :: dha             ! First term A in zero plane displacement height (d/h) (nondimensional)
        real(rk)                   :: dhb             ! Second term B in zero plane displacement height (d/h) (nondimensional)
        real(rk)                   :: qa              ! Coefficient term in Reynolds stress for d/h
        real(rk)                   :: qb              ! Coefficient term in Reynolds stress for d/h
        real(rk)                   :: qc              ! Coefficient term in Reynolds stress for d/h
        real(rk)                   :: qstar           ! Total Reynolds stress term
        real(rk)                   :: fafraczInt_tota ! Numerator term of term B for d/h (Eq. 15 Massman et al.)
        real(rk)                   :: fafraczInt_totb ! Denominator term of term B for d/h (Eq. 15 Massman et al.)
        real(rk)                   :: cstress         ! Suface stress at/above canopy height (nondimensional)
        real(rk)                   :: drag            ! Drag area index (i.e., wind speed attentuation) (nondimensional)
        real(rk)                   :: nrat            ! Ratio of drag/cstress (nondimensional)

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
! Specifically, see Eq. 15 (pg 599) and Eq. 12 (pg 597 of Massman et al.)

        ! Calculate zero-plane displacement height, d/h (Eq. 15 in Massman et al. 2017):
        drag    = CDRAG*PAI
        ustrmod = UBZREF*(0.38 - (0.38 + (vonk/log(Z0GHC)))*exp(-1.0*(15.0*drag)))
        cstress = (2.0*(ustrmod**2.0))/(UBZREF**2.0)
        nrat   =  drag/cstress
        qc = 0.60
        qb = 2.0/(1.0 - exp(-1.0))
        qa = 4.04*(-1.0*qb)
        qstar = qa + qb*exp(-1.0*qc*nrat)
        dha =  1.0 - (cosh(qstar*nrat*FCLAI(1))/cosh(qstar*nrat))
        fafraczInt_tota = IntegrateTrapezoid( ZHC,(cosh(qstar*nrat*FCLAI)*ZHC) )
        fafraczInt_totb = IntegrateTrapezoid( ZHC, cosh(qstar*nrat*FCLAI) )
        dhb = fafraczInt_tota/fafraczInt_totb

        ! zero-plane displacement height
        d_h = dha * dhb
        ! Calculate surface (soil+veg) roughness length, zo/h (Eq. 16 in Massman et al. 2017):
        zo_h  = LAMDARS * (1.0 - d_h) * exp (-vonk*sqrt(2.0/cstress))

    END SUBROUTINE CANOPY_ZPD

end module canopy_zpd_mod
