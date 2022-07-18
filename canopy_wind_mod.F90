module canopy_wind_mod

  implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CANOPY_WIND( HCM, ZK, FAFRACK, UBZREF, Z0GHCM, &
                              CDRAG, PAI, CANBOT_OUT, CANTOP_OUT, CANWIND )

!-----------------------------------------------------------------------

! Description:
!     computes mean wind speed for given height (z) below the canopy top.

! Preconditions:
!     in-canopy height, above-canopy/reference wind speed plant distribution functions

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone canopy wind model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      use canopy_const_mod, ONLY: rk, vonk    !constants for canopy models

! Arguments:
!     IN/OUT
      REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
      REAL(RK),    INTENT( IN )  :: ZK              ! Below canopy height, z (m)
      REAL(RK),    INTENT( IN )  :: FAFRACK         ! Fractional (z) shapes of the
                                                    ! plant surface distribution (nondimensional)
      REAL(RK),    INTENT( IN )  :: UBZREF          ! Mean wind speed at zref-height of canopy top (m/s)
      REAL,        INTENT( IN )  :: Z0GHCM          ! Ratio of ground roughness length to canopy top height (nondimensional)
      REAL(RK),    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
      REAL(RK),    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
      REAL(RK),    INTENT( OUT ) :: CANBOT_OUT      ! Canopy bottom wind reduction factor = canbot (nondimensional)
      REAL(RK),    INTENT( OUT ) :: CANTOP_OUT      ! Canopy top wind reduction factor = cantop    (nondimensional)
      REAL(RK),    INTENT( OUT ) :: CANWIND         ! Mean canopy wind speed at current z (m/s)
!     Local variables
      real(rk)                   :: ustrmod         ! Friction Velocity parameterization (m/s)
      real(rk)                   :: z0g             ! Ground roughness length based on z0g/HCCM ratio (m)
      real(rk)                   :: zkhcm           ! Current zk/hcm ratio (nondimensional)
      real(rk)                   :: cstress         ! Suface stress at/above canopy height (nondimensional)
      real(rk)                   :: drag            ! Drag area index (i.e., wind speed attentuation) (nondimensional)
      real(rk)                   :: nrat            ! Ratio of drag/cstress (nondimensional)
      real(rk)                   :: canbot          ! Logarithmic wind speed that is dominant near the ground (nondimensional)
      real(rk)                   :: cantop          ! Hyperbolic cosine wind speed that is dominant near the top of canopy (nondimensional)

!Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
!(2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
    zkhcm = ZK/HCM
    z0g = Z0GHCM*HCM

    !Nondimensional canopy wind speed term that dominates near the ground:
    if (ZK >= z0g .and. ZK <= HCM) then
      canbot = log(zkhcm/Z0GHCM)/log(1.0_rk/Z0GHCM)
    else if (ZK >= 0 .and. ZK <= z0g) then
      canbot = 0.0  !No-slip condition at surface (u=0 at z=0)
    else
      canbot = 1  ! should be zk > hcm
    end if

    CANBOT_OUT = canbot
    !Nondimensional canopy wind speed term that dominates near the top of the canopy:
    !Assume the drag area distribution over depth of canopy can be approx. p1=0 (no shelter factor) and d1=0
    !(no drag coefficient relation to wind speed) -- thus no intergration then required in Eq. (4) of Massman et al.
    drag    = CDRAG*PAI
    ustrmod = UBZREF*(0.38_rk - (0.38_rk + (vonk/log(Z0GHCM)))*exp(-1.0_rk*(15.0_rk*drag)))
    cstress = (2.0_rk*(ustrmod**2.0_rk))/(UBZREF**2.0_rk)
    nrat   =  drag/cstress
    cantop = cosh(nrat*FAFRACK)/cosh(nrat)
    CANTOP_out = cantop

    if (ZK <= HCM) then
      CANWIND = UBZREF*canbot*cantop
    else
      CANWIND = UBZREF
    end if

    END SUBROUTINE CANOPY_WIND

end module canopy_wind_mod
