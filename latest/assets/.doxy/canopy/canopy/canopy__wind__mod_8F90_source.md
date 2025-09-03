

# File canopy\_wind\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_wind\_mod.F90**](canopy__wind__mod_8F90.md)

[Go to the documentation of this file](canopy__wind__mod_8F90.md)


```Fortran



module canopy_wind_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_wind_most( HCM, ZK, FAFRACK, UBZREF, Z0GHC, &
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
        real(rk)                   :: ustrmod

        real(rk)                   :: z0g

        real(rk)                   :: zkhcm
        real(rk)                   :: cstress

        real(rk)                   :: drag

        real(rk)                   :: nrat

        real(rk)                   :: canbot

        real(rk)                   :: cantop

        real(rk)                   :: zpd

        real(rk)                   :: z0m

        real(rk)                   :: uc

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
        zkhcm = zk/hcm
        z0g = z0ghc*hcm

! Nondimensional canopy wind speed term that dominates near the top of the canopy:
! Assume the drag area distribution over depth of canopy can be approx. p1=0 (no shelter factor) and d1=0
! (no drag coefficient relation to wind speed) -- thus no integration then required in Eq. (4) of Massman et al.
        drag    = cdrag*pai

        zpd = d_h*hcm  !zero-plane displacement height (not scaled to HCM)
        z0m = zo_h*hcm !aerodynamic roughness length (not scaled to HCM)

!First adjust reference wind down to canopy top to get u* from Massman closure equations (necessary) or
!from unified RSL approach.  This is needed to get accurate canopy stress terms for in-canopy winds later...

        if((hcm-zpd) <= 0.) then
            write(*,*) "critical problem: hcan <= zpd"
            call exit(1)
        end if


        if (href > z0m) then ! input wind speed reference height is > roughness length
            uc = ubzref*log(lambdars*(hcm-zpd+z0m)/z0m)/log(href/z0m)  !MOST From NoahMP (M. Barlarge) with user RSL influence term (LAMBDARS)
        else                 ! reference height is <= roughness length--at canopy top (used for observation comparison)
            uc = ubzref
        end if

        !Some checks on Uc calculation
        if (uc > ubzref) then !reference height too small and close to roughness length
            uc = ubzref
        end if

        if (uc <= 0.0) then
            write(*,*)  'Uc cannot be <= 0 at  ',uc , ' in canopy_wind calc...exiting'
            call exit(1)
        end if

        !Calculate U* from Massman 1997 (https://doi.org/10.1023/A:1000234813011)
        ustrmod = uc*(0.38_rk - (0.38_rk + (vonk/log(z0ghc)))*exp(-1.0_rk*(15.0_rk*drag)))

!Calculate in-canopy parameters that dominate near top region of canopy (Massman et al., 2017)
        cstress = (2.0_rk*(ustrmod**2.0_rk))/(uc**2.0_rk)
        nrat   =  drag/cstress
        cantop = cosh(nrat*fafrack)/cosh(nrat)
        cantop_out = cantop

!Calculate in-canopy parameters that dominate near the ground:
        if (zk >= z0g .and. zk <= hcm) then
            canbot = log(zkhcm/z0ghc)/log(1.0_rk/z0ghc)
        else if (zk >= 0 .and. zk <= z0g) then
            canbot = 0.0  !No-slip condition at surface (u=0 at z=0)
        else
            canbot = 1  ! should be zk > hcm
        end if
        canbot_out = canbot

!Calculate the canopy wind variable both above and below canopy top

        if (zk <= hcm) then       !at or below canopy top --> Massman in-canopy profile
            canwind = uc*canbot*cantop
        else                      !above canopy top       --> MOST or RSL profile
            if (uc < ubzref) then !reference height is not small compared to z0m
                canwind = ubzref*log(lambdars*(zk-zpd+z0m)/z0m)/log(href/z0m)
            else                  !cannot calcualate above canopy wind, set constant to UBZREF
                canwind = ubzref
            end if
        end if !Above or below canopy

    END SUBROUTINE canopy_wind_most

end module canopy_wind_mod
```


