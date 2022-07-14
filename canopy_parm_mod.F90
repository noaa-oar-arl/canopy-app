module canopy_parm_mod

  implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CANOPY_PARM( LAT, LON, VTYPE, FCH, FFRAC, LAI, &
                              FIRETYPE, FLAMEH, CDRAG, &
                              PAI, ZCANMAX, SIGMAU, SIGMA1 )

!-----------------------------------------------------------------------
 
! Description:  
!     determines the canopy, fire, and foliage shape distribution parameters
 
! Preconditions:
!     lat, lon, and vegetation type
 
! Subroutines and Functions Called:
 
! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial canopy parameter model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Arguments:
      INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
!     IN/OUT
      REAL(RK),    INTENT( IN )  :: LAT             ! Grid cell latitude (degrees)
      REAL(RK),    INTENT( IN )  :: LON             ! Grid cell longitude (degrees)
      REAL(RK),    INTENT( IN )  :: FCH             ! Grid cell canopy height (m)
      REAL(RK),    INTENT( IN )  :: VTYPE           ! Grid cell dominant vegetation type
      REAL(RK),    INTENT( IN )  :: FFRAC           ! Grid cell forest fraction
      REAL(RK),    INTENT( IN )  :: LAI             ! Grid cell leaf area index

      INTEGER,     INTENT( OUT ) :: FIRETYPE        ! 1 = Above Canopy Fire; 0 = Below Canopy Fire; -1 No Canopy
      REAL(RK),    INTENT( OUT ) :: FLAMEH          ! Flame Height (m)
      REAL(RK),    INTENT( OUT ) :: CDRAG           ! Drag coefficient (nondimensional)
      REAL(RK),    INTENT( OUT ) :: PAI             ! Plant/foliage area index (nondimensional)
      REAL(RK),    INTENT( OUT ) :: ZCANMAX         ! Height of maximum foliage area density (z/h) (nondimensional)
      REAL(RK),    INTENT( OUT ) :: SIGMAU          ! Standard deviation of shape function above zcanmax (z/h) 
      REAL(RK),    INTENT( OUT ) :: SIGMA1          ! Standard deviation of shape function below zcanmax (z/h)

      !Local variables
      !Parameters
         real(rk),    parameter    ::    pi=3.14159    !Value for PI

!need vegtype to Massman forest mapping

!initializing
FIRETYPE=-1
FLAMEH=0.0_rk
CDRAG=0.0_rk
PAI=0.0_rk
ZCANMAX=0.0_rk
SIGMAU=0.0_rk
SIGMA1=0.0_rk

!only populate parameers if contiguous canopy conditions are satistified in each grid cell
if (FCH .gt. 0.5 .and. FFRAC .gt. 0.5 .and. LAI .ge. 0.1) then ! set canopy parameters

        !need better vegtype and lat/lon regional mapping to Massman et al. forest types 
        if (VTYPE .ge. 1 .and. VTYPE .le. 2) then !VIIRS Cat 1-2/Evergreen Needleleaf & Broadleaf 
                                                  !--> Use average Massman Aspen+Spruce+Pine Forest
         FIRETYPE=0
         FLAMEH=2.0_rk
         CDRAG=(0.20_rk + 0.25_rk + 0.20_rk + 0.20_rk + 0.20_rk)/5.0_rk
         PAI=(5.73_rk + 3.28_rk + 2.41_rk + 2.14_rk + 3.78_rk)/5.0_rk  !will use Massman calculation (Eq. 19) of PAI later...
         ZCANMAX=(0.60_rk + 0.36_rk + 0.60_rk + 0.58_rk + 0.60_rk)/5.0_rk
         SIGMAU=(0.38_rk + 0.60_rk + 0.30_rk + 0.20_rk + 0.10_rk)/5.0_rk
         SIGMA1=(0.16_rk + 0.20_rk + 0.10_rk + 0.20_rk + 0.27_rk)/5.0_rk
        end if

        if (VTYPE .ge. 3 .and. VTYPE .le. 5) then !VIIRS Cat 3-5/Deciduous Needleleaf, Broadleaf, Mixed Forests 
                                                  !--> Use Massman Hardwood Forest
         FIRETYPE=0
         FLAMEH=2.0_rk
         CDRAG=0.15_rk
         PAI=4.93_rk !will use Massman calculation (Eq. 19) of PAI later...
         ZCANMAX=0.84_rk
         SIGMAU=0.13_rk
         SIGMA1=0.30_rk
        end if

        if (VTYPE .ge. 6 .and. VTYPE .le. 10 .or. VTYPE .le. 12 ) then !VIIRS Cat 6-10 or 12/Shrubs, Croplands, and Grasses 
                                                     !--> Average of Massman Corn + Rice ) 
         FIRETYPE=1
         FLAMEH=2.0_rk
         CDRAG=(0.30_rk + 0.30_rk)/2.0
         PAI=(2.94_rk + 3.10_rk)/2.0 !will use Massman calculation (Eq. 19) of PAI later...
         ZCANMAX=(0.94_rk + 0.62_rk)/2.0
         SIGMAU=(0.03_rk + 0.50_rk)/2.0
         SIGMA1=(0.60_rk + 0.45_rk)/2.0
        end if


end if
    
    END SUBROUTINE CANOPY_PARM

end module canopy_parm_mod
