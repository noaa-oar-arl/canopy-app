program canopy_driver
!
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
    use canopy_const_mod, ONLY: rk                 !canopy constants
    use canopy_utils_mod, ONLY: IntegrateTrapezoid !utilities for canopy models
    use canopy_parm_mod  !main canopy parameters
    use canopy_files_mod !main canopy input files
    use canopy_wind_mod  !main canopy wind model
    use canopy_waf_mod   !main Wind Adjustment Factor (WAF) model

    implicit none

! !....this block defines geographic domain of inputs (read from user namelist)
    integer        ::    nlat        !length of x coordinate
    integer        ::    nlon        !length of y coordinate
    integer        ::    canlays     !Number of total above and below canopy layers
    real           ::    canres      !Real value of canopy vertical resolution (m)
    real           ::    href        !Reference Height above canopy @ 10 m  (m)
    real           ::    flameh      !Flame Height (m)
    logical        ::    ifcanwind   !logical canopy wind/WAF option (default = .FALSE.)
    integer        ::    pai_opt     !integer for PAI values used or calculated (default = 0)
    real(rk)       ::    pai_set     !real value for PAI set values used (default = 4.0)

! !....this block gives assumed constant parameters for in-canopy conditions (read from user namelist)
    real           ::    z0ghcm   ! ratio of ground roughness length to canopy top height
    !  (default case: Approx. currently does not account for understory
    !   variability)
    real           ::   lamdars     ! Influence function associated with roughness sublayer (nondimensional)

! !....this block gives input canopy height and above reference conditions that should be passed (assume winds at href)
    real(rk)       ::    lat             !latitude  (degrees)
    real(rk)       ::    lon             !longitude (degrees)
    real(rk)       ::    hcm             !Input Canopy Height (m)
    real(rk)       ::    ubzref          !Input above canopy/reference 10-m model wind speed (m/s)
    real(rk)       ::    cluref          !Input canopy clumping index
    real(rk)       ::    lairef          !Input leaf area index
    integer        ::    vtyperef        !Input vegetation type (VIIRS)
    real(rk)       ::    ffracref        !Input forest fraction of grid cell
    real(rk)       ::    ustref          !Input friction velocity
    real(rk)       ::    cszref          !Input cosine of zenith angle
    real(rk)       ::    z0ref           !Input total/surface roughness length
    real(rk)       ::    molref          !Input Monin-Obukhov Length

! !....this block gives vegetion-type canopy dependent  parameters based on Katul et al. (2004)
    integer        ::    firetype         !1 = Above Canopy Fire; 0 = Below Canopy Fire
    real(rk)       ::    cdrag         !Drag coefficient (nondimensional)
    real(rk)       ::    pai           !Plant/foliage area index (nondimensional)
    real(rk)       ::    zcanmax       !Height of maximum foliage area density (z/h) (nondimensional)
    real(rk)       ::    sigmau        !Standard deviation of shape function above zcanmax (z/h)
    real(rk)       ::    sigma1        !Standard deviation of shape function below zcanmax (z/h)

!Local variables
    integer i,i0,loc
    real(rk), allocatable :: zkcm       ( : )  ! in-canopy heights (m)
    real(rk), allocatable :: ztothc     ( : )  ! z/h
    real(rk), allocatable :: fainc      ( : )  ! incremental foliage shape function
    real(rk), allocatable :: fafracz    ( : )  ! incremental fractional foliage shape function
    real(rk), allocatable :: fafraczInt ( : )  ! integral of incremental fractional foliage shape function
    real(rk), allocatable :: canBOT     ( : )  ! Canopy bottom wind reduction factors (nondimensional)
    real(rk), allocatable :: canTOP     ( : )  ! Canopy top wind reduction factors (nondimensional)
    real(rk)              :: fatot             ! integral of total fractional foliage shape function
    real(rk), allocatable :: canWIND    ( :, : )  ! final mean canopy wind speeds (m/s)

    integer  ::    cansublays           ! number of sub-canopy layers
    integer  ::    canmidpoint          ! indice of the sub-canopy midpoint
    integer  ::    flamelays            ! number of flame layers
    integer  ::    midflamepoint        ! indice of the mid-flame point
    real(rk), allocatable :: waf        (:)  ! Calculated Wind Adjustment Factor

!     Test Generic 1D canopy data that should represent virtualized canopy
    TYPE :: profile_type
        integer     :: canlay       !profile layer for model
        real(rk)    :: zk           !profile heights for model (m)
    end TYPE profile_type

    type(profile_type), allocatable :: profile( : )

!     Test Generic 2D met/sfc input variables that should be passed to canopy calculations
    TYPE :: variable_type
        real(rk)   :: lat          !latitude of cell/point
        real(rk)   :: lon          !longitude of cell/point
        real(rk)   :: fh           !forest/canopy height
        real(rk)   :: ws           !wind speed at reference height above canopy (10 m)
        real(rk)   :: clu          !clumping index
        real(rk)   :: lai          !leaf area index
        integer    :: vtype        !vegetation type
        real(rk)   :: ffrac        !forest fraction
        real(rk)   :: ust          !friction velocity (u*)
        real(rk)   :: csz          !cosine of solar zenith angle
        real(rk)   :: z0           !surface roughness length
        real(rk)   :: mol          !Monin-Obukhov length
    end TYPE variable_type

    type(variable_type), allocatable :: variables( : )

!-------------------------------------------------------------------------------
! Read user options from namelist.
!-------------------------------------------------------------------------------

    call  canopy_readnml(nlat,nlon,canlays,canres,href,z0ghcm,lamdars, &
        flameh, ifcanwind, pai_opt, pai_set)
    if (ifcanwind) then
        write(*,*)  'Canopy wind/WAF option selected'
    else
        write(*,*)  'No option(s) selected'
    end if

!-------------------------------------------------------------------------------
! Allocate necessary variables.
!-------------------------------------------------------------------------------

    if(.not.allocated(zkcm)) allocate(zkcm(canlays))
    if(.not.allocated(ztothc)) allocate(ztothc(canlays))
    if(.not.allocated(fainc)) allocate(fainc(canlays))
    if(.not.allocated(fafracz)) allocate(fafracz(canlays))
    if(.not.allocated(fafraczInt)) allocate(fafraczInt(canlays))
    if(.not.allocated(canBOT)) allocate(canBOT(canlays))
    if(.not.allocated(canTOP)) allocate(canTOP(canlays))
    if(.not.allocated(canWIND)) allocate(canWIND(canlays,nlat*nlon))
    if(.not.allocated(waf)) allocate(waf(nlat*nlon))
    if(.not.allocated(profile)) allocate(profile(canlays))
    if(.not.allocated(variables)) allocate(variables(nlat*nlon))

! ... read canopy profile data that should be passed
    open(9,  file=file_prof(1),  status='old')
    i0 = 0
    read(9,*,iostat=i0)  ! skip headline
    do i=1, canlays
        read(9, *) profile(i)
    end do
    close(9)

! ... read met/sfc input variables
    open(8,  file=file_vars(1),  status='old')
    i0 = 0
    read(8,*,iostat=i0)  ! skip headline
    do loc=1, nlat*nlon
        read(8, *) variables(loc)
    end do
    close(8)

    do loc=1, nlat*nlon
        lat      = variables(loc)%lat
        lon      = variables(loc)%lon
        hcm      = variables(loc)%fh
        ubzref   = variables(loc)%ws
        cluref   = variables(loc)%clu
        lairef   = variables(loc)%lai
        vtyperef = variables(loc)%vtype
        ffracref = variables(loc)%ffrac
        ustref   = variables(loc)%ust
        cszref   = variables(loc)%csz
        z0ref    = variables(loc)%z0
        molref   = variables(loc)%mol

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

        call canopy_parm(vtyperef, hcm, ffracref, lairef, &
            pai_opt, pai_set, firetype, cdrag, &
            pai, zcanmax, sigmau, sigma1)

! ... initialize canopy model and integrate to get fractional plant area distribution functions
        zkcm   = profile%zk
        ztothc = zkcm/hcm
        cansublays  = floor(hcm/canres)
        canmidpoint = cansublays/2
        flamelays    = floor(flameh/canres)
        midflamepoint   = flamelays/2

! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
        fainc = 0.0_rk  ! initialize
        fatot = 0.0_rk  ! initialize
        do i=1, canlays
            if (ztothc(i) >= zcanmax .and. ztothc(i) <= 1.0) then
                fainc(i) = exp((-1.0*((ztothc(i)-zcanmax)**2.0))/sigmau**2.0)
            else if (ztothc(i) >= 0.0 .and. ztothc(i) <= zcanmax) then
                fainc(i) = exp((-1.0*((zcanmax-ztothc(i))**2.0))/sigma1**2.0)
            end if
        end do
        fatot = IntegrateTrapezoid(ztothc,fainc)

! ... calculate plant distribution function
        fafracz    = 0.0_rk  ! initialize
        fafraczInt = 0.0_rk  ! initialize
        do i=1, canlays
            fafracz(i) = fainc(i)/fatot
            fafraczInt(i) = IntegrateTrapezoid(ztothc(1:i),fafracz(1:i))
        end do

        if (ifcanwind) then
! ... calculate in-canopy wind speeds at z
            do i=1, canlays
                call canopy_wind(hcm, zkcm(i), fafraczInt(i), ubzref, &
                    z0ghcm, cdrag, pai, canBOT(i), canTOP(i), canWIND(i, loc))
            end do

! ... calculate wind adjustment factor dependent on canopy winds and fire type

            call canopy_waf(hcm, ztothc(1:cansublays), fafraczInt(1:cansublays), fafraczInt(1), &
                ubzref, z0ghcm, lamdars, cdrag, pai, href, flameh, firetype, &
                canBOT(midflamepoint), canTOP(midflamepoint), waf(loc))
        end if

    end do

    if (ifcanwind) then
        write(*,*)  'Writing canopy wind/WAF output'
! ... save as text files for testing
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
    end if

end program canopy_driver
