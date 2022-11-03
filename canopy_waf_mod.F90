module canopy_waf_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_FLAMEH( FLAMEH_OPT, FLAMEH_SET, DX, MODRES, &
        FRP, MIDFLAMEPOINT, FLAMEH )

!-----------------------------------------------------------------------

! Description:
!     computes flame height and midflamepoint layer needed for WAF

! Preconditions:
!     user flameh option and set, grid resolution dx, vert model res, and FRP

! Subroutines and Functions Called:

! Revision History:
!     Oct 2022 P.C. Campbell: Initial flame height calculation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk         !constants for canopy models
        use canopy_utils_mod, ONLY: CalcFlameH !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )  :: FLAMEH_OPT      ! Integer for flameh values used or calculated (default = 0)
        REAL(RK),    INTENT( IN )  :: FLAMEH_SET      ! User Set Flame Height (m)
        REAL(RK),    INTENT( IN )  :: DX              ! DX cell distances using haversine formula (m)
        REAL(RK),    INTENT( IN )  :: MODRES          ! Canopy model input vertical resolution (m)
        REAL(RK),    INTENT( IN )  :: FRP             ! Model input FRP (MW/grid cell area)
        INTEGER,     INTENT( OUT ) :: MIDFLAMEPOINT   ! Indice of the mid-flame point
        REAL(RK),    INTENT( OUT ) :: FLAMEH          ! Flame Height (m)

!     Local variables

        integer  ::    flamelays                      ! number of flame layers

        if (FLAMEH_OPT .eq. 0) then
            if (DX .le. 0.0) then !single lon point -- reset to user value
                write(*,*)  'FLAMEH_OPT set to FRP calc, but dx <= 0...setting flameh = ', &
                    FLAMEH_SET, ' from namelist'
                FLAMEH = FLAMEH_SET
            else                                    !calculate flameh
                FLAMEH = CalcFlameH(FRP,DX)
            end if
        else if (FLAMEH_OPT .eq. 1) then  !user set value
            FLAMEH = FLAMEH_SET
        else
            write(*,*)  'Wrong FLAMEH_OPT choice of ', FLAMEH_OPT, ' in namelist...exiting'
            call exit(2)
        end if
        if (FLAMEH .lt. MODRES) then !flameh beween first (z=0) and second layer height
            MIDFLAMEPOINT = 2    !do not put at z = 0, but rather in second layer
        else
            flamelays     = floor(FLAMEH/MODRES) + 1
            MIDFLAMEPOINT = max(ceiling(flamelays/2.0),2)
        end if

    END SUBROUTINE CANOPY_FLAMEH

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_WAF( HCM, LAMDARS, HREF, FLAMEH, FIRETYPE, &
        D_H, ZO_H, CANBOTMID, CANTOPMID, WAF )

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
        use canopy_const_mod, ONLY: rk !constants for canopy models
        use canopy_utils_mod           !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: LAMDARS         ! Value representing influence of roughness sublayer (nondimensional)
        REAL(RK),    INTENT( IN )  :: HREF            ! Reference Height (m) above the canopy
        REAL(RK),    INTENT( IN )  :: FLAMEH          ! Flame Height (m) -- Only for Above Canopy Fire
        INTEGER ,    INTENT( IN )  :: FIRETYPE        ! 1 = Above Canopy Fire; 0 = Below Canopy Fire
        REAL(RK),    INTENT( IN )  :: CANBOTMID       ! Mid-flame canopy bottom wind reduction factor (nondimensional)
        REAL(RK),    INTENT( IN )  :: CANTOPMID       ! Mid-flame canopy top wind reduction factor (nondimensional)
        REAL(RK),    INTENT( IN )  :: D_H             ! Zero-plane displacement height, d/h
        REAL(RK),    INTENT( IN )  :: ZO_H            ! Surface (soil+veg) roughness length, zo/h
        REAL(RK),    INTENT( OUT ) :: WAF             ! Wind Adjustment Factor (nondimensional)
!     Local variables
        real(rk)                   :: term1           ! Major Term1 in WAF calculation (Eqs. 17 and 18 Massman et al. 2017)
        real(rk)                   :: term2           ! Major Term2 in WAF calculation (Eqs. 17 and 18 Massman et al. 2017)
        real(rk)                   :: delta           ! Ratio parmeter used in WAF for above-canopy (Eq. 18 Massman et al.)

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354

        if (HREF <= 0.0) then
            write(*,*) "critical problem: HREF <= 0, WAF calculation not accurate and thus is reset to 1"
            waf = 1.0
        else
            if (FIRETYPE == 0) then  !sub-canopy
                term1 = log( LAMDARS * ( (1.0 - D_H)/ZO_H ) )  !numerator
                term2 = log( LAMDARS * ( ( (HREF/HCM) + 1.0 - D_H ) / ZO_H ) )  !denominatory
                waf   = CANBOTMID * CANTOPMID * (term1 / term2)
            else                     !above-canopy
                delta = (1.0 - D_H) / (FLAMEH/HCM)
                term1 = log( LAMDARS * ( ( (FLAMEH/HCM) +  1.0 - D_H ) / ZO_H ) ) - &   !numerator
                    1.0 + (delta*log((1.0/delta) + 1.0))
                term2 = log( LAMDARS * ( ( ( (HREF/HCM) + 1.0 - D_H) )/ ZO_H ) )  !denominator
                waf   = term1 / term2
            end if
        end if

    END SUBROUTINE CANOPY_WAF

end module canopy_waf_mod
