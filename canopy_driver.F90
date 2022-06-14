      program canopy_driver

!  the driver to run the canopy winds and WAF algorithms
!
!  W.J. Massman, J.M. Forthofer, and M.A. Finney. An improved 
!  canopy wind model for predicting wind adjustment factors 
!  and wildland fire behavior. Canadian Journal of Forest Research. 
!  47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354
!
!  History:
!    Prototype: Patrick C. Campbell, 06/2022
!
!-------------------------------------------------------------
       use canopy_utils     !utilities for canopy models
       use canopy_wind_mod  !main canopy wind model
!      use canopy_waf_mod   !main canopy wind adjustment factor (WAF) model

      implicit none
        INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
! !....this block gives constant above canopy met conditions that should be passed
        real(rk),    parameter    ::    ubzref=100.0       !Above canopy/reference wind speed (cm/s)

! !....this block gives assumed constant parameters for in-canopy conditions 
        real(rk),    parameter    ::    z0ghccm=0.0025     !! ratio of ground roughness length to canopy top height
                                                       !  (default case: Approx. currently does not account for understory
                                                       !   variability)

! !....this block gives vegetion-type canopy dependent  parameters based on Katul et al. (2004).
!       ***Need to mke this dynamic, e.g., canopy parameter look-up table for regional model application**
!       ***Question:  Can we use GEDI information for canopy heights??***
!       Example:  Hardwood Forest Type (Massman et al. and Katul et al.)
!        integer, parameter        ::    canlays=50         !Number of total canopy layers 
!        real(rk),    parameter    ::    hccm=2250.0        !Canopy Height (m)
!        real(rk),    parameter    ::    cdrag=0.15         !Drag coefficient (nondimensional)
!        real(rk),    parameter    ::    pai=4.93           !Plant/foliage area index (nondimensional)
!        real(rk),    parameter    ::    zcanmax=0.84       !Height of maximum foliage area density (z/h) (nondimensional)
!        real(rk),    parameter    ::    sigmau=0.13        !Standard deviation of shape function above zcanmax (z/h)
!        real(rk),    parameter    ::    sigma1=0.30        !Standard deviation of shape function below zcanmax (z/h)
!       Example:  Aspen Forest Type (Massman et al. and Katul et al.)
        integer, parameter        ::    canlays=50         !Number of total canopy layers
        real(rk),    parameter    ::    hccm=1000.0        !Canopy Height (m)
        real(rk),    parameter    ::    cdrag=0.20         !Drag coefficient (nondimensional)
        real(rk),    parameter    ::    pai=3.28           !Plant/foliage area index (nondimensional)
        real(rk),    parameter    ::    zcanmax=0.36       !Height of maximum foliage area density (z/h) (nondimensional)
        real(rk),    parameter    ::    sigmau=0.60        !Standard deviation of shape function above zcanmax (z/h)
        real(rk),    parameter    ::    sigma1=0.20        !Standard deviation of shape function below zcanmax (z/h)

        integer i,i0
        real(rk) :: zkcm       ( canlays )  ! in-canopy heights (cm)
        real(rk) :: ztothc     ( canlays )  ! z/h
        real(rk) :: fainc      ( canlays )  ! incremental foliage shape function
        real(rk) :: fafracz    ( canlays )  ! incremental fractional foliage shape function
        real(rk) :: fafraczInt ( canlays )  ! integral of incremental fractional foliage shape function
        real(rk) :: canWIND    ( canlays )  ! final mean canopy wind speeds (cm/s)
        real(rk) :: fatot                   ! integral of total fractional foliage shape function

!     met 3D input profile data that should be passed to canopy calculations
      TYPE :: profile_type
           integer     :: canlay       !profile layer for model
           real(rk)    :: zk           !profile heights for model (m)
           real(rk)    :: dzk          !profile height increments (m)
      end TYPE profile_type

      type(profile_type) :: profile( canlays )

! ... read canopy profile data that should be passed
      open(9,  file='input_profile.txt',  status='old')
      i0 = 0
      read(9,*,iostat=i0)  	! skip headline
      do i=1, canlays
        read(9, *) profile(i)
      end do
! ... initialize canopy model and integrate to get fractional plant area distribution functions
      zkcm   = profile%zk*100.0 !convert to cm
      ztothc = zkcm/hccm 
      
 ! calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
      do i=1, canlays
        if (ztothc(i) >= zcanmax .and. ztothc(i) <= 1.0) then
           fainc(i) = exp((-1.0*((ztothc(i)-zcanmax)**2.0))/sigmau**2.0)
        else if (ztothc(i) >= 0.0 .and. ztothc(i) <= zcanmax) then
           fainc(i) = exp((-1.0*((zcanmax-ztothc(i))**2.0))/sigma1**2.0)
        end if
      end do
      fatot      = IntegrateTrapezoid(ztothc,fainc)

! calculate plant distribution function and in-canopy wind speeds at z
      do i=1, canlays
        fafracz(i) = fainc(i)/fatot
        fafraczInt(i) = IntegrateTrapezoid(ztothc(1:i),fafracz(1:i))
        call canopy_wind(hccm, zkcm(i), fafraczInt(i), ubzref, &
                  z0ghccm, cdrag, pai, canWIND(i))
      end do

      write(*,*)  'Below and Above Canopy Wind Speed (cm/s):' , canWIND
      write(*,*)  'Average Wind Speed (cm/s)', SUM(canWIND)/canlays

    end program canopy_driver
