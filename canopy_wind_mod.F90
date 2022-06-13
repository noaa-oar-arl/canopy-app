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
!     IN/OUT
      REAL,    INTENT( IN )  :: HCCM            ! Height of canopy top (cm)
      REAL,    INTENT( IN )  :: ZK              ! Below canopy height, z (cm)
      REAL,    INTENT( IN )  :: FAFRACK         ! Fractional (z) shape of the 
                                                ! plant surface distribution (nondimensional)
      REAL,    INTENT( IN )  :: UBZREF          ! Mean wind speed at zref-height of canopy top (cm/s)
      REAL,    INTENT( IN )  :: Z0GHCCM         ! Ratio of ground roughness length to canopy top height (nondimensional)
      REAL,    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
      REAL,    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
      REAL,    INTENT( OUT ) :: CANWIND         ! Mean canopy wind speed at current z (cm/s)
!     Local calculations
      real(kind=dp)          :: ustrmod         ! Friction Velocity parameterization (cm/s)
      real(kind=dp)          :: z0g             ! Ground roughness length based on z0g/HCCM ratio (cm)
      real(kind=dp)          :: zkhccm          ! Current zk/hccm ratio (nondimensional)
      real(kind=dp)          :: cstress         ! Suface stress at/above canopy height (nondimensional)
      real(kind=dp)          :: drag            ! Drag area index (i.e., wind speed attentuation) (nondimensional)
      real(kind=dp)          :: nrat            ! Ratio of drag/cstress (nondimensional)
      real(kind=dp)          :: canbot          ! Logarithmic wind speed that is dominant near the ground (nondimensional)
      real(kind=dp)          :: cantop          ! Hyperbolic cosine wind speed that is dominant near the top of canopy (nondimensional)

!Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
!(2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
    zkhccm=KZ/HCCM
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

   if (ZK <= HCCM) then
      CANWIND=UBZREF*canbot*cantop
    else
      CANWIND=UBZREF
    end if

end module canopy_wind_mod
