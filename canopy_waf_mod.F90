module canopy_waf_mod

  implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CANOPY_WAF( HCM, ZTOTHC, FAFRACK, FAFRACK0, UBZREF, Z0GHCM, &
                             LAMDARS, CDRAG, PAI, HREF, FLAMEH, FIRETYPE, &
                             CANBOTMID, CANTOPMID, WAF )

!-----------------------------------------------------------------------
 
! Description:  
!     computes Wind Adjustment Factor for fire spread for either sub- or above-canopy fires.
 
! Preconditions:
!     in-canopy height, firetype and height, mean mid-canopy wind speed, and plant distribution functions
 
! Subroutines and Functions Called:
 
! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone canopy wind model 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      use canopy_utils_mod     !utilities for canopy models

! Arguments:
      INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
!     IN/OUT
      REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
      REAL(RK),    INTENT( IN )  :: ZTOTHC  ( : )   ! SUB-CANOPY Total z/h layers (nondimensional)
      REAL(RK),    INTENT( IN )  :: FAFRACK ( : )   ! SUB-CANOPY Fractional (z) shapes of the
                                                    ! plant surface distribution at z=0 (nondimensional)
      REAL(RK),    INTENT( IN )  :: FAFRACK0        ! Fractional (z) shapes of the 
                                                    ! plant surface distribution at z=0 (nondimensional)
      REAL(RK),    INTENT( IN )  :: UBZREF          ! Mean wind speed at zref-height of canopy top (m/s)
      REAL(RK),    INTENT( IN )  :: Z0GHCM          ! Ratio of ground roughness length to canopy top height (nondimensional)
      REAL(RK),    INTENT( IN )  :: LAMDARS         ! Influence function associated with roughness sublayer (nondimensional)
      REAL(RK),    INTENT( IN )  :: HREF            ! Reference Height (m) above the canopy
      REAL(RK),    INTENT( IN )  :: FLAMEH          ! Flame Height (m) -- Only for Above Canopy Fire
      INTEGER ,    INTENT( IN )  :: FIRETYPE        ! 1 = Above Canopy Fire; 0 = Below Canopy Fire
      REAL(RK),    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
      REAL(RK),    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
      REAL(RK),    INTENT( IN )  :: CANBOTMID       ! Mid-flame canopy bottom wind reduction factor (nondimensional)
      REAL(RK),    INTENT( IN )  :: CANTOPMID       ! Mid-flame canopy top wind reduction factor (nondimensional)
      REAL(RK),    INTENT( OUT ) :: WAF             ! Wind Adjustment Factor (nondimensional)
!     Local variables
      ! real(rk)                   :: zkhcm           ! Current zk/hcm ratio (nondimensional)
      real(rk)                   :: ustrmod         ! Friction Velocity parameterization (m/s)
      real(rk)                   :: dha             ! First term A in zero plane displacement height (d/h) (nondimensional)
      real(rk)                   :: dhb             ! Second term B in zero plane displacement height (d/h) (nondimensional)
      real(rk)                   :: d_h             ! Zero displacement height, d/h (i.e., Product of dha*dhb)
      real(rk)                   :: qa              ! Coefficient term in Reynolds stress for d/h
      real(rk)                   :: qb              ! Coefficient term in Reynolds stress for d/h
      real(rk)                   :: qc              ! Coefficient term in Reynolds stress for d/h
      real(rk)                   :: qstar           ! Total Reynolds stress term
      real(rk)                   :: fafraczInt_tota ! Numerator term of term B for d/h (Eq. 15 Massman et al.)
      real(rk)                   :: fafraczInt_totb ! Denominator term of term B for d/h (Eq. 15 Massman et al.)
      real(rk)                   :: cstress         ! Suface stress at/above canopy height (nondimensional)
      real(rk)                   :: drag            ! Drag area index (i.e., wind speed attentuation) (nondimensional)
      real(rk)                   :: nrat            ! Ratio of drag/cstress (nondimensional)
      real(rk)                   :: zo_h            ! Surface (soil+veg) roughness length, zo/h
      real(rk)                   :: term1           ! Major Term1 in WAF calculation (Eqs. 17 and 18 Massman et al. 2017)
      real(rk)                   :: term2           ! Major Term2 in WAF calculation (Eqs. 17 and 18 Massman et al. 2017)
      real(rk)                   :: delta           ! Ratio parmeter used in WAF for above-canopy (Eq. 18 Massman et al.)

!Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
!(2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354

   !Calculate zero plane displacement height, d/h (Eq. 15 in Massman et al. 2017):
   drag    = CDRAG*PAI
   ustrmod = UBZREF*(0.38 - (0.38 + (0.40/log(Z0GHCM)))*exp(-1.0*(15.0*drag)))
   cstress = (2.0*(ustrmod**2.0))/(UBZREF**2.0)
   nrat   =  drag/cstress
   qc = 0.60
   qb = 2.0/(1.0 - exp(-1.0))
   qa = 4.04*(-1.0*qb)
   qstar = qa + qb*exp(-1.0*qc*nrat)
   dha =  1.0 - (cosh(qstar*nrat*FAFRACK0)/cosh(qstar*nrat))

   fafraczInt_tota = IntegrateTrapezoid( ZTOTHC,(cosh(qstar*nrat*FAFRACK)*ZTOTHC) )
   fafraczInt_totb = IntegrateTrapezoid( ZTOTHC, cosh(qstar*nrat*FAFRACK) )
   dhb = fafraczInt_tota/fafraczInt_totb 

   ! zero plane displacement height
   d_h = dha * dhb
   !Calculate surface (soil+veg) roughness length, zo/h (Eq. 16 in Massman et al. 2017):  
   zo_h  = LAMDARS * (1.0 - d_h) * exp (-0.4*sqrt(2.0/cstress))

   !Calculate WAF dependent on fire type (sub- or above-canopy) (Eqs. 17 and 18 of Massman et al. 2017)
   if (FIRETYPE == 0) then  !sub-canopy
    ! write(*,*)  '------Sub-Canopy Fire Type------'
     term1 = log( LAMDARS * ( (1.0 - d_h)/zo_h ) )  !numerator  
     term2 = log( LAMDARS * ( ( (HREF/HCM) + 1.0 - d_h ) / zo_h ) )  !denominatory
     waf   = CANBOTMID * CANTOPMID * (term1 / term2)
   else                     !above-canopy
    ! write(*,*)  '-----Above-Canopy Fire Type-----'
    delta = (1.0 - d_h) / (FLAMEH/HCM)
    term1 = log( LAMDARS * ( ( (FLAMEH/HCM) +  1.0 - d_h ) / zo_h ) ) - &   !numerator
            ( 1.0 + (delta*log((1.0/delta) + 1.0)) )
    term2 = log( LAMDARS * ( ( ( (HREF/HCM) + 1.0 - d_h) )/ zo_h ) )  !denominator
    waf   = term1 / term2
   end if

    END SUBROUTINE CANOPY_WAF

end module canopy_waf_mod
