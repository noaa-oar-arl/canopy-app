module canopy_wind_mod

  implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CANOPY_WIND( HCCM, ZK, FAFRACK, UBZREF, Z0GHCCM, &
                              CDRAG, PAI, CANWIND )

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
! Arguments:
      INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
!     IN/OUT
      REAL(RK),    INTENT( IN )  :: HCCM            ! Height of canopy top (cm)
      REAL(RK),    INTENT( IN )  :: ZK              ! Below canopy height, z (cm)
      REAL(RK),    INTENT( IN )  :: FAFRACK         ! Fractional (z) shapes of the 
                                                ! plant surface distribution (nondimensional)
      REAL(RK),    INTENT( IN )  :: UBZREF          ! Mean wind speed at zref-height of canopy top (cm/s)
      REAL(RK),    INTENT( IN )  :: Z0GHCCM         ! Ratio of ground roughness length to canopy top height (nondimensional)
      REAL(RK),    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
      REAL(RK),    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
      REAL(RK),    INTENT( OUT ) :: CANWIND         ! Mean canopy wind speed at current z (cm/s)
!     Local variables
      real(rk)                   :: ustrmod         ! Friction Velocity parameterization (cm/s)
      real(rk)                   :: z0g             ! Ground roughness length based on z0g/HCCM ratio (cm)
      real(rk)                   :: zkhccm          ! Current zk/hccm ratio (nondimensional)
      real(rk)                   :: cstress         ! Suface stress at/above canopy height (nondimensional)
      real(rk)                   :: drag            ! Drag area index (i.e., wind speed attentuation) (nondimensional)
      real(rk)                   :: nrat            ! Ratio of drag/cstress (nondimensional)
      real(rk)                   :: canbot          ! Logarithmic wind speed that is dominant near the ground (nondimensional)
      real(rk)                   :: cantop          ! Hyperbolic cosine wind speed that is dominant near the top of canopy (nondimensional)

!Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
!(2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
    zkhccm=ZK/HCCM
    z0g = Z0GHCCM*HCCM

   !Nondimensional canopy wind speed term that dominates near the ground:
    if (ZK >= z0g .and. ZK <= HCCM) then
    canbot = log(zkhccm/Z0GHCCM)/log(1.0/Z0GHCCM)
    else if (ZK >= 0.0 .and. ZK <= z0g) then
    canbot = 0.0  !No-slip condition at surface (u=0 at z=0)
    end if

   !Nondimensional canopy wind speed term that dominates near the top of the canopy:
   !Assume the drag area distribution over depth of canopy can be approx. p1=0 (no shelter factor) and d1=0
   !(no drag coefficient relation to wind speed) -- thus no intergration then required in Eq. (4) of Massman et al.
   drag    = CDRAG*PAI
   ustrmod = UBZREF*(0.38 - (0.38 + (0.40/log(Z0GHCCM)))*exp(-1.0*(15.0*drag)))
   cstress = (2.0*(ustrmod**2.0))/(UBZREF**2.0)
   nrat   =  drag/cstress
   cantop = cosh(nrat*FAFRACK)/cosh(nrat)
   print*,cantop
   if (ZK <= HCCM) then
      CANWIND=UBZREF*canbot*cantop
    else
      CANWIND=UBZREF
    end if
    
    END SUBROUTINE CANOPY_WIND

end module canopy_wind_mod
