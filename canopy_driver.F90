program canopy_driver
!
!  The main driver to run the canopy applications
!
!  History:
!    Prototype: Patrick C. Campbell, 06/2022
!
!-------------------------------------------------------------
    use canopy_const_mod, ONLY: rk  !canopy constants
    use canopy_calcdx_mod  !main grid dx calculation
    use canopy_parm_mod    !main canopy parameters
    use canopy_foliage_mod !main canopy foliage distribution
    use canopy_files_mod   !main canopy input files
    use canopy_zpd_mod     !main displacement height model
    use canopy_wind_mod    !main canopy wind model
    use canopy_flameh_mod  !main flame height model
    use canopy_waf_mod     !main Wind Adjustment Factor (WAF) model
    use canopy_eddyx_mod   !main canopy eddy diffusivities
    use canopy_phot_mod    !main canopy photolysis attenuation

    implicit none

! !....this block defines geographic domain of inputs (read from user namelist)
    integer        ::    nlat        !length of x coordinate
    integer        ::    nlon        !length of y coordinate
    integer        ::    modlays     !Number of total above and below canopy model layers
    real(rk)       ::    modres      !Real value of model above and below canopy vertical resolution (m)
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
    !  (default case: Approx. currently does not account for understory variability)
    real(rk)       ::   lamdars      !Value representing influence of roughness sublayer (nondimensional)

! !....this block gives input canopy height and above reference conditions that should be passed (assume winds at href)
    real(rk)       ::    hcmref           !Input Canopy Height (m)
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
    real(rk) ::    flameh                         ! flame Height (m)
    integer  ::    cansublays                     ! number of sub-canopy layers
    integer  ::    midflamepoint                  ! indice of the mid-flame point
!allocatable
    real(rk), allocatable :: zk         ( : )     ! in-canopy heights (m)
    real(rk), allocatable :: zhc        ( : )     ! z/h
    real(rk), allocatable :: fafraczInt ( : )     ! integral of incremental fractional foliage shape function
    real(rk), allocatable :: canBOT     ( : )     ! Canopy bottom wind reduction factors
    real(rk), allocatable :: canTOP     ( : )     ! Canopy top wind reduction factors
    real(rk), allocatable :: canWIND    ( :, : )  ! canopy wind speeds (m/s)
    real(rk), allocatable :: dx         ( : )     ! Model grid cell distance/resolution (m)
    real(rk), allocatable :: waf        ( : )     ! Calculated Wind Adjustment Factor
    real(rk), allocatable :: Kz         ( :, : )  ! Eddy Diffusivities (m2/s)
    real(rk), allocatable :: rjcf       ( :, : )  ! Photolysis Attenuation Correction Factors

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

    call  canopy_readnml(nlat,nlon,modlays,modres,href,z0ghc,lamdars, &
        flameh_opt, flameh_set, ifcanwind, ifcaneddy, ifcanphot,   &
        pai_opt, pai_set, lu_opt, dx_opt, dx_set, lai_thresh, &
        frt_thresh, fch_thresh, rsl_opt)


! ... TODO:  Move these input related things and specific allocations to own subroutine

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

    if(.not.allocated(zk))         allocate(zk(modlays))
    if(.not.allocated(zhc))        allocate(zhc(modlays))
    if(.not.allocated(fafraczInt)) allocate(fafraczInt(modlays))
    if(.not.allocated(canBOT))     allocate(canBOT(modlays))
    if(.not.allocated(canTOP))     allocate(canTOP(modlays))
    if(.not.allocated(canWIND))    allocate(canWIND(modlays,nlat*nlon))
    if(.not.allocated(dx))         allocate(dx(nlat*nlon))
    if(.not.allocated(waf))        allocate(waf(nlat*nlon))
    if(.not.allocated(Kz))         allocate(Kz(modlays,nlat*nlon))
    if(.not.allocated(rjcf))       allocate(rjcf(modlays,nlat*nlon))
    if(.not.allocated(variables))  allocate(variables(nlat*nlon))

! ... TODO:  NL condition for input file/data type read from txt or netcdf 1D or 2D

! ... read met/sfc input variables from text file
    open(8,  file=file_vars(1),  status='old')
    i0 = 0
    read(8,*,iostat=i0)  ! skip headline
    do loc=1, nlat*nlon
        read(8, *) variables(loc)
    end do
    close(8)

! ... TODO:  add option to read met/sfc input variables from 1D ncf file

! ... TODO:  add option to read met/sfc input variables from 2D ncf file

! ... derive canopy model profile heights from user NL input
    zk(1) = 0.0_rk
    do i=2, modlays
        zk(i)   = zk(i-1) + modres
    end do

! ... initialize ***grid cell_only*** dependent variables
    dx   = 0.0_rk
    waf  = 1.0_rk

! ... get dx cell distances using haversine formula

! ... TODO:  Set nlat = nlat_user (if txt) or nlat = nlat_file (if ncf)
! ... TODO:  Set nlon = nlat_user (if txt) or nlon = nlon_file (if ncf)

    call canopy_calcdx(dx_opt, dx_set, nlat, nlon, variables%lat, &
        variables%lon, dx)

! ... Main loop through model grid cells
    do loc=1, nlat*nlon
        hcmref   = variables(loc)%fh
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

! ... get scaled canopy model profile and sub-canopy layers
        zhc         = zk/hcmref
        cansublays  = floor(hcmref/modres)

! ... initialize ***canopy profile*** dependent variables
        fafraczInt        = 0.0_rk
        canTOP            = 1.0_rk
        canBOT            = 1.0_rk

! ... initialize ***grid cell_and_canopy*** profile dependent variables
        canWIND(:,loc)    = ubzref    !initialize to above canopy wind
        Kz(:,loc)         = -999.0_rk !initialize to missing
        rjcf(:,loc)       = 1.0_rk      !initialize to above canopy rjcf=1

! ... check for model vegetation types
        if (lu_opt .eq. 0 ) then !VIIRS
            if (vtyperef .le. 10 .or. vtyperef .eq. 12) then !VIIRS types

! ... check for contiguous canopy conditions at each model grid cell
                if (hcmref .gt. fch_thresh .and. ffracref .gt. frt_thresh &
                    .and. lairef .gt. lai_thresh) then

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

                    call canopy_parm(vtyperef, hcmref, ffracref, lairef, &
                        pai_opt, pai_set, lu_opt, firetype, cdrag, &
                        pai, zcanmax, sigmau, sigma1)

! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z

                    call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                        fafraczInt)

! ... calculate zero-plane displacement height/hc and surface (soil+veg) roughness lengths/hc

                    call canopy_zpd(zhc(1:cansublays), fafraczInt(1:cansublays), &
                        ubzref, z0ghc, lamdars, cdrag, pai, d_h, zo_h)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                    if (ifcanwind) then
                        do i=1, modlays
                            call canopy_wind(hcmref, zk(i), fafraczInt(i), ubzref, &
                                z0ghc, cdrag, pai, href, d_h, zo_h, molref, &
                                rsl_opt, canBOT(i), canTOP(i), canWIND(i, loc))
                        end do
! ... determine midflamepoint and flame height from user or FRP calculation
                        call canopy_flameh(flameh_opt, flameh_set, dx(loc), modres, &
                            frpref, midflamepoint, flameh)

                        if (flameh .gt. 0.0) then !only calculate WAF when flameh > 0
                            call canopy_waf(hcmref, lamdars, href, flameh, firetype, d_h, zo_h, &
                                canBOT(midflamepoint), canTOP(midflamepoint), waf(loc))
                        end if !flameh_opt
                    end if !canwind

! ... user option to calculate in-canopy eddy diffusivities at height z
                    if (ifcaneddy) then
                        do i=1, modlays
                            call canopy_eddyx(hcmref, zk(i), ustref, molref, Kz(i, loc))
                        end do
                    end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                    if (ifcanphot) then
                        call canopy_phot(fafraczInt, &
                            lairef, cluref, cszref, rjcf(:, loc))
                    end if

                end if !Contiguous Canopy
            end if   !Vegetation types
        else
            write(*,*)  'Wrong LU_OPT choice of ', lu_opt, ' in namelist...exiting'
            call exit(2)

        end if       !Landuse Options
    end do

! ... TODO: add separate subroutine for handling output files (text or ncf based on user NL input)

    if (ifcanwind) then
        write(*,*)  'Writing canopy wind/WAF output'
! ... save as text files for testing
        open(10, file='output_canopy_wind.txt')
        write(10, '(a30, f6.1, a2)') 'Reference height, h: ', href, 'm'
        write(10, '(a30, i6)') 'Number of model layers: ', modlays
        write(10, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'WS (m/s)'
        do loc=1, nlat*nlon
            do i=1, modlays
                write(10, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                    zk(i), canWIND(i, loc)
            end do
        end do

        open(11, file='output_waf.txt')
        write(11, '(a30, f6.1)') 'Reference height, h: ', href, 'm'
        write(11, '(a8, a9, a19, a11)') 'Lat', 'Lon', 'Canopy height (m)', 'WAF'
        do loc=1, nlat*nlon
            write(11, '(f8.2, f9.2, f19.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, hcmref, waf(loc)
        end do
    end if

    if (ifcaneddy) then
        write(*,*)  'Writing canopy eddy diffusivity output'
! ... save as text files for testing
        open(12, file='output_eddy_Kz.txt')
        write(12, '(a30, f6.1, a2)') 'Reference height, h: ', href, 'm'
        write(12, '(a30, i6)') 'Number of model layers: ', modlays
        write(12, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'Kz'
        do loc=1, nlat*nlon
            do i=1, modlays
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
        write(13, '(a30, i6)') 'Number of model layers: ', modlays
        write(13, '(a8, a9, a12, a15)') 'Lat', 'Lon', 'Height (m)', 'rjcf'
        do loc=1, nlat*nlon
            do i=1, modlays
                write(13, '(f8.2, f9.2, f12.2, es15.7)')  variables(loc)%lat, variables(loc)%lon, &
                    zk(i), rjcf(i, loc)
            end do
        end do
    end if

! ... TODO:  Add deallocating subroutine/submodule for specific variables allocated

    if(allocated(zk))         deallocate(zk)
    if(allocated(zhc))        deallocate(zhc)
    if(allocated(fafraczInt)) deallocate(fafraczInt)
    if(allocated(canBOT))     deallocate(canBOT)
    if(allocated(canTOP))     deallocate(canTOP)
    if(allocated(canWIND))    deallocate(canWIND)
    if(allocated(dx))         deallocate(dx)
    if(allocated(waf))        deallocate(waf)
    if(allocated(Kz))         deallocate(Kz)
    if(allocated(rjcf))       deallocate(rjcf)
    if(allocated(variables))  deallocate(variables)

end program canopy_driver
