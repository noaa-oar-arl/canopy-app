module canopy_wind_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_WIND_MOST( HCM, ZK, FAFRACK, UBZREF, Z0GHC, &
        CDRAG, PAI, HREF, D_H, ZO_H, LAMBDARS, &
        CANBOT_OUT, CANTOP_OUT, CANWIND )

!-----------------------------------------------------------------------

! Description:
!     computes mean wind speed for given height (z) below the canopy top.

! Preconditions:
!     in-canopy height, at canopy wind speed, plant distribution functions

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     below the canopy and use of MOST above canopy.
!     Jun 2022 P.C. Campbell: Initial standalone canopy wind model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk, vonk  !constants for canopy models

! Arguments:
!       IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM              ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: ZK               ! Above/Below canopy height, z (m)
        REAL(RK),    INTENT( IN )  :: FAFRACK          ! Fractional (z) shapes of the
        ! plant surface distribution (nondimensional)
        REAL(RK),    INTENT( IN )  :: UBZREF           ! Mean wind speed at reference height (ZREF) (m/s)
        REAL(RK),    INTENT( IN )  :: Z0GHC            ! Ratio of ground roughness length to canopy top height (nondimensional)
        REAL(RK),    INTENT( IN )  :: LAMBDARS         ! Value representing influence of roughness sublayer (nondimensional)
        REAL(RK),    INTENT( IN )  :: CDRAG            ! Drag coefficient (nondimensional)
        REAL(RK),    INTENT( IN )  :: PAI              ! Total plant/foliage area index (nondimensional)
        REAL(RK),    INTENT( IN )  :: HREF             ! Reference Height above the canopy (ZREF = HCM + HREF; m)
        REAL(RK),    INTENT( IN )  :: D_H              ! Zero plane displacement heights (nondimensional)
        REAL(RK),    INTENT( IN )  :: ZO_H             ! Surface (soil+veg) roughness lengths (nondimensional)
        REAL(RK),    INTENT( OUT ) :: CANBOT_OUT       ! Canopy bottom wind reduction factor = canbot (nondimensional)
        REAL(RK),    INTENT( OUT ) :: CANTOP_OUT       ! Canopy top wind reduction factor = cantop    (nondimensional)
        REAL(RK),    INTENT( OUT ) :: CANWIND          ! Mean canopy wind speed at current z (m/s)
!       Local variables
        real(rk)                   :: ustrmod          ! Friction Velocity parameterization based on Massman 2017 (m/s)
        real(rk)                   :: z0g              ! Ground roughness length based on z0g/HCCM ratio (m)
        real(rk)                   :: zkhcm            ! Current zk/hcm ratio (nondimensional)
        real(rk)                   :: cstress          ! Surface stress at/above canopy height (nondimensional)
        real(rk)                   :: drag             ! Drag area index (i.e., wind speed attenuation) (nondimensional)
        real(rk)                   :: nrat             ! Ratio of drag/cstress (nondimensional)
        real(rk)                   :: canbot           ! Logarithmic wind speed that is dominant near the ground (nondimensional)
        real(rk)                   :: cantop           ! Hyperbolic cosine wind speed that is dominant near the top of canopy (nondimensional)
        real(rk)                   :: zpd              ! Zero plane displacement heights MOST  (m)
        real(rk)                   :: z0m              ! Surface (soil+veg) roughness lengths with Massman LAMBDARS (m)
        real(rk)                   :: uc               ! Wind directly at canopy top (m/s)

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
        zkhcm = ZK/HCM
        z0g = Z0GHC*HCM

! Nondimensional canopy wind speed term that dominates near the top of the canopy:
! Assume the drag area distribution over depth of canopy can be approx. p1=0 (no shelter factor) and d1=0
! (no drag coefficient relation to wind speed) -- thus no integration then required in Eq. (4) of Massman et al.
        drag    = CDRAG*PAI

        zpd = D_H*HCM  !zero-plane displacement height (not scaled to HCM)
        z0m = ZO_H*HCM !aerodynamic roughness length (not scaled to HCM)

!First adjust reference wind down to canopy top to get u* from Massman closure equations (necessary) or
!from unified RSL approach.  This is needed to get accurate canopy stress terms for in-canopy winds later...

        if((HCM-zpd) <= 0.) then
            write(*,*) "critical problem: hcan <= zpd"
            call exit(1)
        end if

        if (HREF > z0m) then ! input wind speed reference height is > roughness length
            uc = (vonk/(0.38_rk - (0.38_rk + (vonk/log(Z0GHC)))*exp(-1.0_rk*(15.0_rk*drag)))) * &
                (UBZREF/log((LAMBDARS*(HCM-zpd))/z0m))  !MOST Log Profile from combining M17 Eqs. 10 14
        else                 ! reference height is <= roughness length--at canopy top (used for observation comparison)
            uc = UBZREF
        end if

        !Some checks on Uc calculation
        if (uc > UBZREF) then !reference height too small and close to roughness length
            uc = UBZREF
        end if

        if (uc <= 0.0) then
            write(*,*)  'Uc cannot be <= 0 at  ',uc , ' in canopy_wind calc...exiting'
            call exit(1)
        end if

        !Calculate U* from Massman 1997 (https://doi.org/10.1023/A:1000234813011)
        ustrmod = uc*(0.38_rk - (0.38_rk + (vonk/log(Z0GHC)))*exp(-1.0_rk*(15.0_rk*drag)))

!Calculate in-canopy parameters that dominate near top region of canopy (Massman et al., 2017)
        cstress = (2.0_rk*(ustrmod**2.0_rk))/(uc**2.0_rk)
        nrat   =  drag/cstress
        cantop = cosh(nrat*FAFRACK)/cosh(nrat)
        CANTOP_out = cantop

!Calculate in-canopy parameters that dominate near the ground:
        if (ZK >= z0g .and. ZK <= HCM) then
            canbot = log(zkhcm/Z0GHC)/log(1.0_rk/Z0GHC)
        else if (ZK >= 0 .and. ZK <= z0g) then
            canbot = 0.0  !No-slip condition at surface (u=0 at z=0)
        else
            canbot = 1  ! should be zk > hcm
        end if
        CANBOT_OUT = canbot

!Calculate the canopy wind variable both above and below canopy top

        if (ZK <= HCM) then       !at or below canopy top --> Massman in-canopy profile
            CANWIND = uc*canbot*cantop
        else                      !above canopy top       --> MOST or RSL profile
            if (uc < UBZREF) then !reference height is not small compared to z0m
                CANWIND = UBZREF*log(LAMBDARS*(ZK-zpd+z0m)/z0m)/log(HREF/z0m)
            else                  !cannot calcualate above canopy wind, set constant to UBZREF
                CANWIND = UBZREF
            end if
        end if !Above or below canopy

    END SUBROUTINE CANOPY_WIND_MOST

end module canopy_wind_mod
