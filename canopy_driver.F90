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
    use canopy_utils_mod, ONLY: IntegrateTrapezoid,CalcFlameH !utilities for canopy models
    use canopy_parm_mod    !main canopy parameters
    use canopy_files_mod   !main canopy input files
    use canopy_zpd_mod     !main displacement height model
    use canopy_wind_mod    !main canopy wind model
    use canopy_waf_mod     !main Wind Adjustment Factor (WAF) model
    use canopy_eddyx_mod   !main canopy eddy diffusivities
    use canopy_phot_mod   !main canopy photolysis attenuation

    implicit none

! !....this block defines geographic domain of inputs (read from user namelist)
    integer        ::    nlat        !length of x coordinate
    integer        ::    nlon        !length of y coordinate
    integer        ::    canlays     !Number of total above and below canopy layers
    real(rk)       ::    canres      !Real value of canopy vertical resolution (m)
    real(rk)       ::    href        !Reference Height above canopy @ 10 m  (m)
    logical        ::    ifcanwind   !logical canopy wind/WAF option (default = .FALSE.)
    logical        ::    ifcaneddy   !logical canopy eddy Kz option (default = .FALSE.)
    logical        ::    ifcanphot   !logical canopy photolsyis atten option (default = .FALSE.)
    integer        ::    pai_opt     !integer for PAI values used or calculated (default = 0)
    real(rk)       ::    pai_set     !real value for PAI set values used (default = 4.0)
    integer        ::    lu_opt      !integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
    integer        ::    flameh_opt  !Integer for flameh values used or calculated (default = 0)
    real(rk)       ::    flameh_set  !User Set Flame Height (m)
    integer        ::    dx_opt      !Integer for dx resolution values used or calculated (default = 0)
    real(rk)       ::    dx_set      !User Set Grid Cell Resolution (m)
    real(rk)       ::    lai_thresh  !User set grid cell LAI threshold to apply canopy conditions (m2/m2)
    real(rk)       ::    frt_thresh  !User set grid cell forest fraction threshold to apply canopy conditions ()
    real(rk)       ::    fch_thresh  !User set grid cell canopy height threshold to apply canopy conditions (m)
    integer        ::    rsl_opt     !RSL option used in model from Rosenzweig et al. 2021 (default = 0, off)


! !....this block gives assumed constant parameters for in-canopy conditions (read from user namelist)
    real(rk)       ::   z0ghc   ! ratio of ground roughness length to canopy top height
    !  (default case: Approx. currently does not account for understory
    !   variability)
    real(rk)       ::   lamdars     ! Influence function associated with roughness sublayer (nondimensional)

! !....this block gives input canopy height and above reference conditions that should be passed (assume winds at href)
    real(rk)       ::    lat             !latitude  (degrees)
    real(rk)       ::    lon             !longitude (degrees)
    real(rk)       ::    dlon            !longitude increment (degrees)
    real(rk)       ::    hcm             !Input Canopy Height (m)
    real(rk)       ::    dx              !Model grid cell resolution (west-east) (m)
    real(rk)       ::    ubzref          !Input above canopy/reference 10-m model wind speed (m/s)
    real(rk)       ::    cluref          !Input canopy clumping index
    real(rk)       ::    lairef          !Input leaf area index
    integer        ::    vtyperef        !Input vegetation type (VIIRS)
    real(rk)       ::    ffracref        !Input forest fraction of grid cell
    real(rk)       ::    ustref          !Input friction velocity
    real(rk)       ::    cszref          !Input cosine of zenith angle
    real(rk)       ::    z0ref           !Input total/surface roughness length
    real(rk)       ::    molref          !Input Monin-Obukhov Length
    real(rk)       ::    frpref          !Input fire radiative power

! !....this block gives vegetion-type canopy dependent  parameters based on Katul et al. (2004)
    integer        ::    firetype      !1 = Above Canopy Fire; 0 = Below Canopy Fire
    real(rk)       ::    cdrag         !Drag coefficient (nondimensional)
    real(rk)       ::    pai           !Plant/foliage area index (nondimensional)
    real(rk)       ::    zcanmax       !Height of maximum foliage area density (z/h) (nondimensional)
    real(rk)       ::    sigmau        !Standard deviation of shape function above zcanmax (z/h)
    real(rk)       ::    sigma1        !Standard deviation of shape function below zcanmax (z/h)
    real(rk)       ::    d_h           !Zero plane displacement heights (z/h)
    real(rk)       ::    zo_h          !Surface (soil+veg) roughness lengths (z/h)
!Local variables
    integer i,i0,loc
    real(rk), allocatable :: zk         ( : )     ! in-canopy heights (m)
    real(rk), allocatable :: zhc     ( : )     ! z/h
    real(rk), allocatable :: fainc      ( : )     ! incremental foliage shape function
    real(rk), allocatable :: fafracz    ( : )     ! incremental fractional foliage shape function
    real(rk), allocatable :: fafraczInt ( : )     ! integral of incremental fractional foliage shape function
    real(rk), allocatable :: canBOT     ( : )     ! Canopy bottom wind reduction factors
    real(rk), allocatable :: canTOP     ( : )     ! Canopy top wind reduction factors
    real(rk)              :: fatot                ! integral of total fractional foliage shape function
    real(rk), allocatable :: canWIND    ( :, : )  ! canopy wind speeds (m/s)
    real(rk), allocatable :: Kz         ( :, : )  ! Eddy Diffusivities
    real(rk), allocatable :: RJCF       ( :, : )  ! Photolysis Attenuation Coefficients

    real(rk) ::    flameh               ! flame Height (m)
    integer  ::    cansublays           ! number of sub-canopy layers
    integer  ::    flamelays            ! number of flame layers
    integer  ::    midflamepoint        ! indice of the mid-flame point
    real(rk), allocatable :: waf        (:)  ! Calculated Wind Adjustment Factor

! Generic 2D met/sfc input variables that should be passed to canopy calculations
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
        real(rk)   :: frp          !fire radiative power (MW per grid cell)
    end TYPE variable_type

    type(variable_type), allocatable :: variables( : )

!-------------------------------------------------------------------------------
! Read user options from namelist.
!-------------------------------------------------------------------------------

    call  canopy_readnml(nlat,nlon,canlays,canres,href,z0ghc,lamdars, &
        flameh_opt, flameh_set, ifcanwind, ifcaneddy, ifcanphot,   &
        pai_opt, pai_set, lu_opt, dx_opt, dx_set, lai_thresh, &
        frt_thresh, fch_thresh, rsl_opt)

    if (ifcanwind) then
        write(*,*)  'Canopy wind/WAF option selected'
    end if
    if (ifcaneddy) then
        write(*,*)  'Canopy eddy Kz option selected'
    end if
    if (ifcanphot) then
        write(*,*)  'Canopy photolysis attenuation option selected'
    end if
    if (.not. any([ifcanwind, ifcaneddy, ifcanphot])) then
        write(*,*)  'No option(s) selected...see namelist'
        call exit(0)
    end if

!-------------------------------------------------------------------------------
! Allocate necessary variables.
!-------------------------------------------------------------------------------

    if(.not.allocated(zk)) allocate(zk(canlays))
    if(.not.allocated(zhc)) allocate(zhc(canlays))
    if(.not.allocated(fainc)) allocate(fainc(canlays))
    if(.not.allocated(fafracz)) allocate(fafracz(canlays))
    if(.not.allocated(fafraczInt)) allocate(fafraczInt(canlays))
    if(.not.allocated(canBOT)) allocate(canBOT(canlays))
    if(.not.allocated(canTOP)) allocate(canTOP(canlays))
    if(.not.allocated(canWIND)) allocate(canWIND(canlays,nlat*nlon))
    if(.not.allocated(waf)) allocate(waf(nlat*nlon))
    if(.not.allocated(Kz)) allocate(Kz(canlays,nlat*nlon))
    if(.not.allocated(RJCF)) allocate(RJCF(canlays,nlat*nlon))
    if(.not.allocated(variables)) allocate(variables(nlat*nlon))

! ... read met/sfc input variables
    open(8,  file=file_vars(1),  status='old')
    i0 = 0
    read(8,*,iostat=i0)  ! skip headline
    do loc=1, nlat*nlon
        read(8, *) variables(loc)
    end do
    close(8)

! ... if multiple grid points, calculate grid cell east-west distance
    if (dx_opt .eq. 0) then !user set to calculate dx grid cell distance from grid lons
        if (nlon .gt. 1) then !convert lon grid points to distance (m)
            dlon = abs(variables(nlat*nlon)%lon - variables((nlat*nlon)-1)%lon)
            dx   = dlon*111000.0_rk
        else                  !single grid cell/point, use namelist defined dx resolution (m) for cell
            write(*,*)  'DX_OPT set to calc, but nlon <= 1...setting dx = ', &
                dx_set, ' from namelist'
            dlon = 0.0_rk
            dx = dx_set
        end if
    else ! user set dx_set from namelist
        dx = dx_set
    end if

! ... derive canopy model profile heights
    zk(1) = 0.0_rk
    do i=2, canlays
        zk(i)   = zk(i-1) + canres
    end do

! ... initialize grid cell dependent variables
    waf  = 1.0_rk

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
        frpref   = variables(loc)%frp

! ... get scaled canopy model profile and layers
        zhc      = zk/hcm
        cansublays  = floor(hcm/canres)

! ... initialize canopy profile dependent variables
        fainc             = 0.0_rk
        fatot             = 0.0_rk
        fafracz           = 0.0_rk
        fafraczInt        = 0.0_rk
        canTOP            = 1.0_rk
        canBOT            = 1.0_rk

! ... initialize grid cell and canopy profile dependent variables
        canWIND(:,loc)    = ubzref    !initialize to above canopy wind
        Kz(:,loc)         = -999.0_rk !initialize to missing
        RJCF(:,loc)     = 1.0_rk    !initialize to above canopy RJCF=1

! ... check for model vegetation types
        if (vtyperef .le. 10 .or. vtyperef .eq. 12) then

! ... check for contiguous canopy conditions at each model grid cell
            if (hcm .gt. fch_thresh .and. ffracref .gt. frt_thresh &
                .and. lairef .ge. lai_thresh) then

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

                call canopy_parm(vtyperef, hcm, ffracref, lairef, &
                    pai_opt, pai_set, lu_opt, firetype, cdrag, &
                    pai, zcanmax, sigmau, sigma1)

! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
                do i=1, canlays
                    if (zhc(i) >= zcanmax .and. zhc(i) <= 1.0) then
                        fainc(i) = exp((-1.0*((zhc(i)-zcanmax)**2.0))/sigmau**2.0)
                    else if (zhc(i) >= 0.0 .and. zhc(i) <= zcanmax) then
                        fainc(i) = exp((-1.0*((zcanmax-zhc(i))**2.0))/sigma1**2.0)
                    end if
                end do
                fatot = IntegrateTrapezoid(zhc,fainc)

! ... calculate plant distribution function
                do i=1, canlays
                    fafracz(i) = fainc(i)/fatot
                    fafraczInt(i) = IntegrateTrapezoid(zhc(1:i),fafracz(1:i))
                end do

! ... calculate zero-plane displacement height/hc and surface (soil+veg) roughness lengths/hc

                call canopy_zpd(zhc(1:cansublays), fafraczInt(1:cansublays), &
                    ubzref, z0ghc, lamdars, cdrag, pai, d_h, zo_h)

! ... user option to calculate in-canopy wind speeds at height z

                if (ifcanwind) then

                    do i=1, canlays
                        call canopy_wind(hcm, zk(i), fafraczInt(i), ubzref, &
                            z0ghc, cdrag, pai, href, d_h, zo_h, molref, &
                            rsl_opt, canBOT(i), canTOP(i), canWIND(i, loc))
                    end do

! ... calculate wind adjustment factor dependent on fires (FRP), canopy winds,
!                                                        flameh, and fire type
                    if (flameh_opt .eq. 0) then
                        if (dx .le. 0.0) then !single lon point -- reset to user value
                            write(*,*)  'FLAMEH_OPT set to FRP calc, but dx <= 0...setting flameh = ', &
                                flameh_set, ' from namelist'
                            flameh = flameh_set
                        else                                    !calculate flameh
                            flameh = CalcFlameH(frpref,dx)
                        end if
                    else if (flameh_opt .eq. 1) then  !user set value
                        flameh = flameh_set
                    else
                        write(*,*)  'Wrong FLAMEH_OPT choice of ', flameh_opt, ' in namelist...exiting'
                        call exit(2)
                    end if
                    if (flameh .lt. canres) then !flameh under first layer height
                        midflamepoint = 2 !put in layer above
                    else
                        flamelays     = floor(flameh/canres)
                        midflamepoint = max((flamelays/2),2)
                    end if
                    if (flameh .gt. 0.0) then !only calculate when flameh > 0
                        call canopy_waf(hcm, lamdars, href, flameh, firetype, d_h, zo_h, &
                            canBOT(midflamepoint), canTOP(midflamepoint), waf(loc))
                    end if
                end if

! ... user option to calculate in-canopy eddy diffusivities at height z
                if (ifcaneddy) then
                    do i=1, canlays
                        call canopy_eddyx(hcm, zk(i), ustref, molref, Kz(i, loc))
                    end do
                end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                if (ifcanphot) then
                    call canopy_phot(fafraczInt, &
                        lairef, cluref, cszref, RJCF(:, loc))
                end if

            end if
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
                write(10, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                    zk(i), canWIND(i, loc)
            end do
        end do

        open(11, file='output_waf.txt')
        write(11, '(a30, f6.1)') 'Reference height, h: ', href, 'm'
        write(11, '(a8, a9, a19, a11)') 'Lat', 'Lon', 'Canopy height (m)', 'WAF'
        do loc=1, nlat*nlon
            write(11, '(f8.2, f9.2, f19.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, variables(loc)%fh, waf(loc)
        end do
    end if

    if (ifcaneddy) then
        write(*,*)  'Writing canopy eddy diffusivity output'
! ... save as text files for testing
        open(12, file='output_eddy_Kz.txt')
        write(12, '(a30, f6.1, a2)') 'Reference height, h: ', href, 'm'
        write(12, '(a30, i6)') 'Number of in-canopy layers: ', canlays
        write(12, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'Kz'
        do loc=1, nlat*nlon
            do i=1, canlays
                write(12, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                    zk(i), Kz(i, loc)
            end do
        end do
    end if

    if (ifcanphot) then
        write(*,*)  'Writing canopy photolysis attenuation output'
! ... save as text files for testing
        open(13, file='output_phot.txt')
        write(13, '(a30, f6.1, a2)') 'Reference height, h: ', href, 'm'
        write(13, '(a30, i6)') 'Number of in-canopy layers: ', canlays
        write(13, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'RJCF'
        do loc=1, nlat*nlon
            do i=1, canlays
                write(13, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                    zk(i), RJCF(i, loc)
            end do
        end do
    end if


end program canopy_driver
