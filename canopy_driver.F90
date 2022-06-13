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

! !....this block gives constant above canopy met conditions that should be passed
        real,    parameter    ::    ubzref=100.0       !Above canopy/reference wind speed (cm/s)

! !....this block gives assumed constant parameters for in-canopy conditions 
        real,    parameter    ::    z0ghccm=0.0025     !! ratio of ground roughness length to canopy top height
                                                       !  (default case: Approx. currently does not account for understory
                                                       !   variability)

! !....this block gives vegetion-type canopy dependent  parameters based on Katul et al. (2004).
!       ***Need to mke this dynamic, e.g., canopy parameter look-up table for regional model application**
!       ***Question:  Can we use GEDI information for canopy heights??***
!       Example:  Hardwood Forest Type (Massman et al. and Katul et al.)
        integer, parameter    ::    canlays=35         !Number of total canopy layers 
        real,    parameter    ::    hccm=250.0         !Canopy Height (cm)
        real,    parameter    ::    cdrag=0.15         !Drag coefficient (nondimensional)
        real,    parameter    ::    pai=4.93           !Plant/foliage area index (nondimensional)
        real,    parameter    ::    zcanmax=0.84       !Height of maximum foliage area density (z/h) (nondimensional)
        real,    parameter    ::    sigmau=0.13        !Standard deviation of shape function above zcanmax (z/h)
        real,    parameter    ::    sigma1=0.30        !Standard deviation of shape function below zcanmax (z/h)

        integer i,i0
        real :: canWIND  ( canlays )  ! final mean canopy wind speeds (cm/s)

!     met 3D input profile data that should be passed to canopy calculations
      TYPE :: profile_type
        real    :: zk           !below canopy heights for model (cm)
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

 ! calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
      do i=1, canlays
        ztothc(i) = z(i)/hccm
        if (ztothc(i) >= zmaxrho .and. ztothc(i) <= 1.0) then
           fainc(i) = exp((-1.0*((ztothc(i)-zmaxrho)**2.0))/sigmau**2.0)
        else if (ztothc(i) >= 0.0 .and. ztothc(i) <= zmaxrho) then
           fainc(i) = exp((-1.0*((zmaxrho-ztothc(i))**2.0))/sigma1**2.0)
        end if
      end do
      fatot      = IntegrateTrapezoid(ztothc,fainc)

! calculate plant distribution function and in-canopy wind speeds at z
      do i=1, canlays
        fafracz(i) = fainc(i)/fatot
        fafraczInt(i) = IntegrateTrapezoid(ztothc(1:i),fafracz(1:i))
        
        call canopy_wind(hccm, profile%zk(i), fafraczInt(i), ubzref, &
                  z0ghccm, cdrag, pai, canWIND(i))
      end do

      write(*,*)  'Canopy Wind Speed (cm/s):' , canWIND

    end program canopy_driver
