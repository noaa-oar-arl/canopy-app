      program canopy_driver

!  the driver to run the canopy applications
!
!  Current Applications:
!  1. Canopy winds and Wind Adjustment Factor
!  Citation(s)
!  W.J. Massman, J.M. Forthofer, and M.A. Finney. An improved 
!  canopy wind model for predicting wind adjustment factors 
!  and wildland fire behavior. Canadian Journal of Forest Research. 
!  47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354
!
!  History:
!    Prototype: Patrick C. Campbell, 06/2022
!
!-------------------------------------------------------------
       use canopy_utils_mod !utilities for canopy models
       use canopy_wind_mod  !main canopy wind model
       use canopy_waf_mod   !main Wind Adjustment Factor (WAF) model

      implicit none
        INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
! !....this block defines geographic domain of inputs (CONUS, 25-65N, 50-150W, 0.1 degree resolution)
        integer, parameter        ::    nlat=401        !length of x coordinate
        integer, parameter        ::    nlon=1001       !length of y coordinate
        integer, parameter        ::    canlays=100     !Number of total above and below canopy layers
        real(rk),    parameter    ::    href=10.0       !Reference Height above canopy @ 10 m  (m)

! !....this block gives input canopy height and above reference conditions that should be passed (assume winds at href)
        real(rk)                  ::    hcm             !Input Canopy Height (m)
        real(rk)                  ::    ubzref          !Input above canopy/reference 10-m model wind speed (m/s)

! !....this block gives assumed constant parameters for in-canopy conditions (Unavoidable Uncertainties!!!)
!      Need to move to module canopy_parm_mod for all canopy/vegetation types
        real(rk),    parameter    ::    z0ghcm=0.0025   ! ratio of ground roughness length to canopy top height
                                                        !  (default case: Approx. currently does not account for understory
                                                        !   variability)
        real(rk),    parameter    ::   lamdars=1.25     ! Influence function associated with roughness sublayer (nondimensional)
! !....this block gives vegetion-type canopy dependent  parameters based on Katul et al. (2004)
!      Need to move to module canopy_parm_mod dependent on canopy/vegetation type inputs
!       Example:  Hardwood Forest Type (Massman et al. and Katul et al.)
!        integer, parameter        ::    firetype=0         !1 = Above Canopy Fire; 0 = Below Canopy Fire
!        real(rk),    parameter    ::    flameh=2.0         !Flame Height (m) 
!        real(rk),    parameter    ::    cdrag=0.15         !Drag coefficient (nondimensional)
!        real(rk),    parameter    ::    pai=4.93           !Plant/foliage area index (nondimensional)
!        real(rk),    parameter    ::    zcanmax=0.84       !Height of maximum foliage area density (z/h) (nondimensional)
!        real(rk),    parameter    ::    sigmau=0.13        !Standard deviation of shape function above zcanmax (z/h)
!        real(rk),    parameter    ::    sigma1=0.30        !Standard deviation of shape function below zcanmax (z/h)
!       Example:  Aspen Forest Type (Massman et al. and Katul et al.)
!        integer, parameter        ::    firetype=0         !1 = Above Canopy Fire; 0 = Below Canopy Fire
!        real(rk),    parameter    ::    flameh=2.0         !Flame Height (m) 
!        real(rk),    parameter    ::    cdrag=0.20         !Drag coefficient (nondimensional)
!        real(rk),    parameter    ::    pai=3.28           !Plant/foliage area index (nondimensional)
!        real(rk),    parameter    ::    zcanmax=0.36       !Height of maximum foliage area density (z/h) (nondimensional)
!        real(rk),    parameter    ::    sigmau=0.60        !Standard deviation of shape function above zcanmax (z/h)
!        real(rk),    parameter    ::    sigma1=0.20        !Standard deviation of shape function below zcanmax (z/h)
!       Example:  Corn Crop Type (Massman et al. and Katul et al.)
        integer, parameter        ::    firetype=1         !1 = Above Canopy Fire; 0 = Below Canopy Fire
        real(rk),    parameter    ::    flameh=2.0         !Flame Height (m) 
        real(rk),    parameter    ::    cdrag=0.30         !Drag coefficient (nondimensional)
        real(rk),    parameter    ::    pai=2.94           !Plant/foliage area index (nondimensional)
        real(rk),    parameter    ::    zcanmax=0.94       !Height of maximum foliage area density (z/h) (nondimensional)
        real(rk),    parameter    ::    sigmau=0.03        !Standard deviation of shape function above zcanmax (z/h)
        real(rk),    parameter    ::    sigma1=0.60        !Standard deviation of shape function below zcanmax (z/h)

!Local variables        
        integer i,i0,loc
        real(rk) :: zkcm       ( canlays )  ! in-canopy heights (m)
        real(rk) :: resz       ( canlays )  ! canopy height resolution (m)
        real(rk) :: ztothc     ( canlays )  ! z/h
        real(rk) :: fainc      ( canlays )  ! incremental foliage shape function
        real(rk) :: fafracz    ( canlays )  ! incremental fractional foliage shape function
        real(rk) :: fafraczInt ( canlays )  ! integral of incremental fractional foliage shape function
        real(rk) :: canBOT     ( canlays )  ! Canopy bottom wind reduction factors (nondimensional)
        real(rk) :: canTOP     ( canlays )  ! Canopy top wind reduction factors (nondimensional)
        real(rk) :: fatot                   ! integral of total fractional foliage shape function
        real(rk) :: canWIND    ( canlays, nlat*nlon )  ! final mean canopy wind speeds (m/s)        

        integer  ::    cansublays           ! number of sub-canopy layers
        integer  ::    canmidpoint          ! indice of the sub-canopy midpoint
        integer  ::    flamelays            ! number of flame layers
        integer  ::    midflamepoint        ! indice of the mid-flame point
        real(rk) ::    waf     (nlat*nlon)  ! Calculated Wind Adjustment Factor

!     Test Generic 1D canopy data that should be passed to canopy calculations
      TYPE :: profile_type
           integer     :: canlay       !profile layer for model
           real(rk)    :: zk           !profile heights for model (m)
           real(rk)    :: dzk          !profile height increments (m)
      end TYPE profile_type

      type(profile_type) :: profile( canlays )
      
!     Test Generic 2D met/sfc input variables that should be passed to canopy calculations
      TYPE :: variable_type
           real(rk)    :: lat          !latitude of cell/point
           real(rk)    :: lon          !longitude of cell/point
           real(rk)    :: fh           !forest/canopy height
           real(rk)    :: ws           !wind speed at reference height above canopy (10 m)
      end TYPE variable_type

      type(variable_type) :: variables(nlat*nlon)      


! ... read canopy profile data that should be passed
      open(9,  file='input_profile.txt',  status='old')
      i0 = 0
      read(9,*,iostat=i0)  ! skip headline
      do i=1, canlays
        read(9, *) profile(i)
      end do
! ... read canopy height and reference 10-m wind data
      open(8,  file='input_variables.txt',  status='old')
      i0 = 0
      read(8,*,iostat=i0)  ! skip headline
      do loc=1, nlat*nlon
        read(8, *) variables(loc)
      end do


      do loc=1, nlat*nlon
        hcm    = variables(loc)%fh
        ubzref = variables(loc)%ws

! ... initialize canopy model and integrate to get fractional plant area distribution functions
        zkcm   = profile%zk 
        ztothc = zkcm/hcm 
        resz   = profile%dzk
        cansublays  = floor(hcm/resz(canlays))
        canmidpoint = cansublays/2
        flamelays    = floor(flameh/resz(canlays))
        midflamepoint   = flamelays/2
! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
        fainc = 0  ! initialize
        do i=1, canlays
          if (ztothc(i) >= zcanmax .and. ztothc(i) <= 1.0) then
             fainc(i) = exp((-1.0*((ztothc(i)-zcanmax)**2.0))/sigmau**2.0)
          else if (ztothc(i) >= 0.0 .and. ztothc(i) <= zcanmax) then
             fainc(i) = exp((-1.0*((zcanmax-ztothc(i))**2.0))/sigma1**2.0)
          end if
        end do
        fatot = IntegrateTrapezoid(ztothc,fainc)

! ... calculate plant distribution function and in-canopy wind speeds at z
        do i=1, canlays
          fafracz(i) = fainc(i)/fatot
          fafraczInt(i) = IntegrateTrapezoid(ztothc(1:i),fafracz(1:i))
          call canopy_wind(hcm, zkcm(i), fafraczInt(i), ubzref, &
                    z0ghcm, cdrag, pai, canBOT(i), canTOP(i), canWIND(i, loc))
        end do

!        write(*, '(a10, f7.1, X, a10, f7.1)')  'Latitude:', variables(loc)%lat, 'Longitude:', variables(loc)%lon
!        write(*,*)  'Below and Above (10-m) Canopy Wind Speeds (m/s):'
!        write(*,'(a6, a6, a15)') 'z (m)', 'z/hc', 'WS (m/s)'
!        do i = 1 , canLays
!          write(*,'(f6.2, f6.2, es15.7)') zkcm(i), ztothc(i), canWIND(i, loc)
!        end do
      
        call canopy_waf(hcm, ztothc(1:cansublays), fafraczInt(1:cansublays), fafraczInt(1), & 
                        ubzref, z0ghcm, lamdars, cdrag, pai, href, flameh, firetype, & 
                        canBOT(midflamepoint), canTOP(midflamepoint), waf(loc))
      
!            write(*,*)  'Wind Adjustment Factor:', waf(loc)
      end do


! ... save as text file
      open(10, file='output_canopy_wind.txt')
      write(10, '(a30, f6.1, a2)') 'Reference height, h: ', href, 'm'
      write(10, '(a30, i6)') 'Number of in-canopy layers: ', canlays
      write(10, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'WS (m/s)'
      do loc=1, nlat*nlon
        do i=1, canlays
          write(10, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, profile(i)%zk, canWIND(i, loc)
        end do
      end do

      open(11, file='output_waf.txt')
      write(11, '(a30, f6.1)') 'Reference height, h: ', href, 'm'
      write(11, '(a8, a9, a19, a11)') 'Lat', 'Lon', 'Canopy height (m)', 'WAF'
      do loc=1, nlat*nlon
        write(11, '(f8.2, f9.2, f19.2, f11.7)')  variables(loc)%lat, variables(loc)%lon, variables(loc)%fh, waf(loc)
      end do

    end program canopy_driver
