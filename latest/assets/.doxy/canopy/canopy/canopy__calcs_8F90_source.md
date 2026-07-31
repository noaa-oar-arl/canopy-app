

# File canopy\_calcs.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_calcs.F90**](canopy__calcs_8F90.md)

[Go to the documentation of this file](canopy__calcs_8F90.md)


```Fortran



SUBROUTINE canopy_calcs(nn)

    use canopy_const_mod      !> constants for canopy models
    use canopy_coord_mod      !> main canopy coordinate descriptions
    use canopy_canopts_mod    !> main canopy option descriptions
    use canopy_canmet_mod     !> main canopy met/sfc input descriptions
    use canopy_canvars_mod    !> main canopy variables descriptions
    use canopy_utils_mod      !> main canopy utilities
    use canopy_dxcalc_mod     !> main canopy dx calculation
    use canopy_profile_mod    !> main canopy foliage profile routines
    use canopy_var3din_mod    !> main canopy 3d variable in routines
    use canopy_rad_mod        !> main canopy radiation sunlit/shaded routines
    use canopy_tleaf_mod      !> main canopy leaf temperature sunlit/shaded routines
    use canopy_wind_mod       !> main canopy wind components
    use canopy_fire_mod       !> fire-related canopy calculations
    use canopy_phot_mod       !> photolysis attenuation calculations
    use canopy_eddy_mod       !> eddy diffusivity calculations
    use canopy_bioemi_mod     !> biogenic emission calculations
    use canopy_drydep_mod     !> gas dry deposition calculations
    use canopy_aero_ddep_mod  !> aerosol dry deposition calculations

    IMPLICIT NONE

    INTEGER,     INTENT( IN )       :: nn

    integer i,j,k,loc
    INTEGER  :: int_nlaic
    INTEGER,  save :: int_nlaip
    REAL(rk) :: nlaic, nlaip
    REAL(rk), save :: pastlai
    REAL(rk), save :: currentlai
    REAL(rk) :: tsteplai
    REAL(rk) :: dnewfrac,doldfrac,hnewfrac,holdfrac
    REAL(rk) :: rib,ra
    REAL(rk) :: lat2d(nlon,nlat), lon2d(nlon,nlat)
    REAL(rk) :: lat1d(nlon*nlat), lon1d(nlon*nlat)

    write(*,*)  'Calculating Canopy Parameters'
    write(*,*)  '-------------------------------'

    lat2d = variables_2d%lat
    lon2d = variables_2d%lon
    lat1d = variables%lat
    lon1d = variables%lon

    if (infmt_opt .eq. 0) then !Main input format is 2D NetCDF and output will be 2D NetCDf

        if (ifcanwind .or. ifcanwaf) then !only calculate if canopy wind or WAF option
            call canopy_calcdx_2d(dx_opt, dx_set, nlat, nlon, lat2d, &
                lon2d, dx_2d)
        end if

        if (href_opt .eq. 0 ) then !setting entire array = href_set value from user NL
            variables_2d%href = href_set
        else if (href_opt .eq. 1 ) then !from file array
            variables_2d%href =  variables_2d%href
        else
            write(*,*)  'Wrong HREF_OPT choice of ', href_opt, ' in namelist...exiting'
            call exit(2)
        end if

! ... Main loop through model grid cells for 2D input/output

        do i=1, nlon
            do j=1, nlat

                hcmref         = variables_2d(i,j)%ch
                uref           = variables_2d(i,j)%ugrd10m
                vref           = variables_2d(i,j)%vgrd10m
                cluref         = variables_2d(i,j)%clu
                lairef         = variables_2d(i,j)%lai
                vtyperef       = variables_2d(i,j)%vtype
                canfracref     = variables_2d(i,j)%canfrac
                ustref         = variables_2d(i,j)%fricv
                cszref         = variables_2d(i,j)%csz
                z0ref          = variables_2d(i,j)%sfcr
                molref         = variables_2d(i,j)%mol
                frpref         = variables_2d(i,j)%frp
                hgtref         = variables_2d(i,j)%href
                sotypref       = variables_2d(i,j)%sotyp
                pressfcref     = variables_2d(i,j)%pressfc
                dswrfref       = variables_2d(i,j)%dswrf
                shtflref       = variables_2d(i,j)%shtfl
                tmpsfcref      = variables_2d(i,j)%tmpsfc
                tmp2mref       = variables_2d(i,j)%tmp2m
                spfh2mref      = variables_2d(i,j)%spfh2m
                hpblref        = variables_2d(i,j)%hpbl
                prate_averef   = variables_2d(i,j)%prate_ave
                soilw1ref      = variables_2d(i,j)%soilw1
                soilw2ref      = variables_2d(i,j)%soilw2
                soilw3ref      = variables_2d(i,j)%soilw3
                soilw4ref      = variables_2d(i,j)%soilw4
                wiltref        = variables_2d(i,j)%wilt
                ozone_w126ref  = variables_2d(i,j)%ozone_w126
                soilt1ref      = variables_2d(i,j)%soilt1
                soilt2ref      = variables_2d(i,j)%soilt2
                soilt3ref      = variables_2d(i,j)%soilt3
                soilt4ref      = variables_2d(i,j)%soilt4
                tmp_hyblev1ref = variables_2d(i,j)%tmp_hyblev1
                snowc_averef   = variables_2d(i,j)%snowc_ave
                icec_averef    = variables_2d(i,j)%icec
                vbdsf_averef   = variables_2d(i,j)%vbdsf_ave
                vddsf_averef   = variables_2d(i,j)%vddsf_ave

! ... calculate wind speed from u and v
                ubzref   = sqrt((uref**2.0) + (vref**2.0))

! ... get scaled canopy model profile and sub-canopy layers
                if (hcmref > 0) then
                    zhc = zk / hcmref
                else
                    zhc = 0
                end if
                cansublays  = min(floor(hcmref/modres),modlays)
                if (cansublays .lt. 1) then !case where model resolution >= canopy height
                    cansublays=1            !only one layer allowed
                end if
! ... calculate initial canopy temp/pressure/humidity/density profiles from first guess approximations (i.e., no leaf energy balance)

                do k=1, modlays
                    !For temperature we use 1st model layer and 2-m temperatures to estimate lapse rate, not tmpsfc
                    !which is too extreme in some global locations.
                    tka_3d(i,j,k)=calctemp(zk(k)*100.0_rk, (hyblev1*100.0_rk) - 200.0_rk, &
                        tmp_hyblev1ref-273.15_rk, tmp2mref-273.15_rk) ! temp    [K]
                    pressa_3d(i,j,k)=calcpressure(zk(k)*100.0_rk, 200.0_rk, pressfcref*0.01_rk, &
                        tmp2mref, tmpsfcref)                                                                  ! press   [mb]
                    relhuma_3d(i,j,k)=calcrelhum(tka_3d(i,j,k),pressa_3d(i,j,k),spfh2mref*1000.0_rk)    ! relhum  [%]
                    spechuma_3d(i,j,k)=calcspechum(relhuma_3d(i,j,k),tka_3d(i,j,k),pressa_3d(i,j,k))     ! spechum [g/kg]
                end do

! ... check for valid model vegetation types
                if (lu_opt .eq. 0 .or. lu_opt .eq. 1 ) then !VIIRS or MODIS
                    if (vtyperef .gt. 0 .and. vtyperef .le. 12 &
                        .or. vtyperef .eq. 14 .or. vtyperef .eq. 18 &
                        .or. vtyperef .eq. 19) then

! ... check for can_opt from user namelist

                        if (vtyperef .ge. 1 .and. vtyperef .le. 5 &
                            .or. vtyperef .eq. 18 .or. vtyperef .eq. 19) then !!VIIRS/MODIS forest/tundra canopies
                            if (can_opt .eq. 0) then !use inputs for canopy heights
                                hcmref = hcmref
                                canfracref = canfracref
                                lairef = lairef
                            else if (can_opt .eq. 1) then !user set constant canopy heights for forests/tundras
                                !only used to overide valid vtypes with hcm<=0, canfrac<=0, or lai<=0.
                                if (hcmref .le. 0.0) then
                                    hcmref = can_chset
                                else
                                    hcmref = hcmref
                                endif
                                if (canfracref .le. 0.0) then
                                    canfracref = can_cfset
                                else
                                    canfracref = canfracref
                                endif
                                if (lairef .le. 0.0) then
                                    lairef = can_laiset
                                else
                                    lairef = lairef
                                endif
                                !recalculate
                                zhc         = zk/hcmref
                                cansublays  = min(floor(hcmref/modres),modlays)
                                if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                    cansublays=1            !only one layer allowed
                                end if
                            else
                                write(*,*)  'Wrong CAN_OPT choice of ', can_opt, &
                                    ' in namelist...exiting'
                                call exit(2)
                            end if
                        end if

! ... check for ssg_opt from user namelist
                        if (vtyperef .ge. 6 .and. vtyperef .le. 11) then !VIIRS/MODIS shrubs/savannas/grasses (SSG) type
                            !includes wetlands!
                            if (ssg_opt .eq. 0) then !use inputs for SSG (possibly not measured...)
                                hcmref = hcmref
                                canfracref = canfracref
                                lairef = lairef
                            else if (ssg_opt .eq. 1) then !user set constant shrubs/savannas/grasslands height
                                !only used to overide valid vtypes with hcm<=0, canfrac<=0, or lai<=0.
                                if (hcmref .le. 0.0) then
                                    hcmref = ssg_chset
                                else
                                    hcmref = hcmref
                                endif
                                if (canfracref .le. 0.0) then
                                    canfracref = ssg_cfset
                                else
                                    canfracref = canfracref
                                endif
                                if (lairef .le. 0.0) then
                                    lairef = ssg_laiset
                                else
                                    lairef = lairef
                                endif
                                !recalculate
                                zhc         = zk/hcmref
                                cansublays  = min(floor(hcmref/modres),modlays)
                                if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                    cansublays=1            !only one layer allowed
                                end if
                            else
                                write(*,*)  'Wrong SSG_OPT choice of ', ssg_opt, &
                                    ' in namelist...exiting'
                                call exit(2)
                            end if
                        end if

! ... check for crop_opt from user namelist
                        if (vtyperef .eq. 12 .or. vtyperef .eq. 14) then !VIIRS/MODIS crop types
                            if (crop_opt .eq. 0) then !use inputs for crop (possibly not measured...)
                                hcmref = hcmref
                                canfracref = canfracref
                                lairef = lairef
                            else if (crop_opt .eq. 1) then !user set constant crop height
                                !only used to overide valid vtypes with hcm<=0, canfrac<=0, or lai<=0.
                                if (hcmref .le. 0.0) then
                                    hcmref = crop_chset
                                else
                                    hcmref = hcmref
                                endif
                                if (canfracref .le. 0.0) then
                                    canfracref = crop_cfset
                                else
                                    canfracref = canfracref
                                endif
                                if (lairef .le. 0.0) then
                                    lairef = crop_laiset
                                else
                                    lairef = lairef
                                endif
                                !recalculate
                                zhc         = zk/hcmref
                                cansublays  = min(floor(hcmref/modres),modlays)
                                if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                    cansublays=1            !only one layer allowed
                                end if
                            else
                                write(*,*)  'Wrong CROP_OPT choice of ', crop_opt, &
                                    ' in namelist...exiting'
                                call exit(2)
                            end if
                        end if

! ... check for contiguous canopy conditions at each model grid cell
                        if (hcmref .gt. ch_thresh .and. canfracref .gt. cf_thresh &
                            .and. lairef .gt. lai_thresh) then

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

                            call canopy_parm(vtyperef, hcmref, canfracref, lairef, &
                                pai_opt, pai_set, lu_opt, firetype, cdrag, &
                                pai, zcanmax, sigmau, sigma1)

! ... Choose between prescribed canopy/foliate shape profile or observed GEDI PAVD profile
                            if (pavd_opt .eq. 0) then
! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
                                call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                    fafraczint)
                            else
                                pavdref = variables_3d(i,j,:)%pavd
                                levref  = variables_1d%lev
! ... derive canopy/foliage distribution shape profile from interpolated GEDI PAVD profile - bottom up total in-canopy and fraction at z
                                if (variables_2d(i,j)%lat .gt. (-1.0_rk*pavd_set) .and. &
                                    variables_2d(i,j)%lat .lt. pavd_set) then !use GEDI PAVD
                                    call canopy_pavd2fafrac(zcanmax, sigmau, sigma1, hcmref, zhc, &
                                        pavdref, levref, fafraczint)
                                    !check if there is observed canopy height but no PAVD profile
                                    if (hcmref .gt. 0.0 .and. maxval(fafraczint) .le. 0.0) then !revert to prescribed shape profile
                                        call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                            fafraczint)
                                    end if
                                else !revert back to using prescribed shape profile
                                    call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                        fafraczint)
                                end if
                            end if

! ... calculate leaf area density profile from foliage shape function for output (m2/m3)
                            do k=1, modlays
                                if (zk(k) .gt. 0.0 .and. zk(k) .le. hcmref) then ! above ground level and at/below canopy top
                                    if (k .lt. modlays)  then
                                        lad_3d(i,j,k) = ((fafraczint(k+1) - fafraczint(k))*lairef)/modres
                                    else
                                        lad_3d(i,j,k) = lad_3d(i,j,modlays-1)
                                    end if
                                else
                                    lad_3d(i,j,k) = 0.0_rk
                                end if
                            end do

! ... calculate zero-plane displacement height/hc and surface (soil+veg) roughness lengths/hc
                            call canopy_zpd(zhc(1:cansublays), fafraczint(1:cansublays), &
                                ubzref, z0ghc, lambdars, cdrag, pai, hcmref, hgtref, &
                                z0ref, vtyperef, lu_opt, z0_opt, d_h_2d(i,j), zo_h_2d(i,j))

! ... calculate canopy radiation (sunlit and shade) profile

                            call canopy_fsun_clu( fafraczint, lairef, cluref, cszref, fsun)

! ... calculate canopy leaf temperature (sun/shade) profile

                            call canopy_tleaf_lin(zk, hcmref, tmp2mref, fsun, &
                                tleaf_sun, tleaf_shade, tleaf_ave)

! ... calculate canopy Photosynthetic Photon Flux Density (PPFD) (sun/shade) profile

                            call canopy_ppfd_exp(zk, hcmref, vbdsf_averef, vddsf_averef, lairef, fsun, &
                                ppfd_sun, ppfd_shade, ppfd_ave)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                            if (ifcanwind .or. ifcanwaf) then
                                if (rsl_opt .eq. 0) then
                                    do k=1, modlays
                                        call canopy_wind_most(hcmref, zk(k), fafraczint(k), ubzref, &
                                            z0ghc, cdrag, pai, hgtref, d_h_2d(i,j), zo_h_2d(i,j), &
                                            lambdars, canbot(k), cantop(k), canwind_3d(i,j,k))
                                    end do
                                else
                                    write(*,*) 'wrong RSL_OPT namelist option = ', rsl_opt, 'only option = 0 right now'
                                    call exit(2)
                                end if

! ... determine midflamepoint and flame height from user or FRP calculation
                                call canopy_flameh(flameh_opt, flameh_set, dx_2d(i,j), modres, &
                                    frpref, frp_fac, hcmref, lu_opt, vtyperef, flameh_cal, &
                                    midflamepoint, flameh_2d(i,j))
                                if (firetype .eq. 0) then !forest/sub-canopy firetype
                                    if (flameh_2d(i,j) .gt. 0.0) then !flameh must be > 0
                                        if (flameh_2d(i,j) .le. hcmref) then !only calculate when flameh <= FCH
                                            call canopy_waf(hcmref, lambdars, hgtref, flameh_2d(i,j), &
                                                firetype, d_h_2d(i,j), zo_h_2d(i,j), canbot(midflamepoint), &
                                                cantop(midflamepoint), waf_2d(i,j))
                                        else
                                            write(*,*) 'warning...sub-canopy type fire, but flameh > FCH, setting WAF=1'
                                            waf_2d(i,j) = 1.0_rk
                                        end if
                                    end if
                                else  !grass/crops, above-canopy firetype
                                    if (flameh_2d(i,j) .gt. 0.0) then !flameh still must be > 0
                                        call canopy_waf(hcmref, lambdars, hgtref, flameh_2d(i,j), &
                                            firetype, d_h_2d(i,j), zo_h_2d(i,j), canbot(midflamepoint), &
                                            cantop(midflamepoint), waf_2d(i,j))
                                    end if
                                end if
                                if (waf_2d(i,j) .gt. 1.0_rk) then !Final check of WAF > 1, must be <=1
                                    waf_2d(i,j) = 1.0_rk
                                end if
                            end if

! ... user option to calculate in-canopy eddy diffusivities at height z
                            if (ifcaneddy) then
                                do k=1, modlays
                                    call canopy_eddyx(hcmref, zk(k), ustref, molref, kz_3d(i,j,k))
                                end do
                            end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                            if (ifcanphot) then
                                if (cszref .ge. 0.0_rk) then !only calculate if cell isn't dark
                                    call canopy_phot(fafraczint, &
                                        lairef, cluref, cszref, rjcf_3d(i,j,:))
                                end if
                            end if

!.......user option to calculate in-canopy leafage influence and assigning LAI as per timestep

                            if (leafage_opt .eq. 0) then
                                !!! Check if the lai_tstep is greater than time_intvl
                                if (lai_tstep .ge. time_intvl) then
                                    tsteplai = lai_tstep/86400.0_rk  !convert lai_tstep to time step in # of days needed for bud break in leafage
                                    nlaic = nn*time_intvl/lai_tstep  !calculate number of lai tsteps that have elapsed
                                    int_nlaic = int(nlaic)
                                else
                                    WRITE (*, *) "Error: Input LAI time step cannot be less than model time step...exiting!!!"
                                    CALL exit(1)
                                endif
                                ! Initialize pastlai and currentlai based on current timestep
                                if (nn .eq. 1) then
                                    currentlai = lairef
                                    pastlai = currentlai
                                else if (int_nlaic .gt. int_nlaip) then !only update with each new lai tstep
                                    pastlai = currentlai
                                    currentlai = lairef
                                endif
                                nlaip = nlaic !set value to compare in next model time step
                                int_nlaip = int(nlaip)
                            end if !leafage_opt = 0 end

!.......user option to calculate historical leaf temperature and PAR for past 24-hours and 240-hours rolling average per timestep
                            !Initialize

                            if (hist_opt .eq. 0) then      !Use instantaneous only
                                ppfd_sun24_3d(i,j,:)     = ppfd_sun
                                ppfd_shade24_3d(i,j,:)   = ppfd_shade
                                tleaf_sun24_3d(i,j,:)    = tleaf_sun
                                tleaf_shade24_3d(i,j,:)  = tleaf_shade
                                tleaf_ave24_3d(i,j,:)    = tleaf_ave
                                ppfd_sun240_3d(i,j,:)    = ppfd_sun
                                ppfd_shade240_3d(i,j,:)  = ppfd_shade
                                tleaf_sun240_3d(i,j,:)   = tleaf_sun
                                tleaf_shade240_3d(i,j,:) = tleaf_shade
                                tleaf_ave240_3d(i,j,:)   = tleaf_ave
                                !for AQ, temp, and wind stress factors in biogenics, use instantaneous only (caution!)
                                daily_maxt2m_2d(i,j)     = tmp2mref
                                daily_mint2m_2d(i,j)     = tmp2mref
                                daily_maxws10m_2d(i,j)   = ubzref
                            else if (hist_opt .eq. 1) then  !Try for historical average values
                                ! Calculate weights for running means of historic variables
                                ! DNEWFRAC and DOLDFRAC are the weights given to the current
                                ! and existing value, respectively, when updating running means
                                ! over the last X days. HNEWFRAC and HOLDFRAC are the same but
                                ! for the 24H means.
                                dnewfrac = time_intvl / ( tau_days * 24.0_rk * 3600.0_rk )
                                doldfrac = 1.0_rk - dnewfrac
                                hnewfrac = time_intvl / ( tau_hours * 3600.0_rk )
                                holdfrac = 1.0_rk - hnewfrac

                                if (nn .le. 24) then
                                    !Track times for moving time window (running) average below
                                    ppfd_sun24_tmp_3d(nn,i,j,:)     = ppfd_sun
                                    ppfd_shade24_tmp_3d(nn,i,j,:)   = ppfd_shade
                                    tleaf_sun24_tmp_3d(nn,i,j,:)    = tleaf_sun
                                    tleaf_shade24_tmp_3d(nn,i,j,:)  = tleaf_shade
                                    tleaf_ave24_tmp_3d(nn,i,j,:)    = tleaf_ave
                                    ppfd_sun240_tmp_3d(nn,i,j,:)    = ppfd_sun
                                    ppfd_shade240_tmp_3d(nn,i,j,:)  = ppfd_shade
                                    tleaf_sun240_tmp_3d(nn,i,j,:)   = tleaf_sun
                                    tleaf_shade240_tmp_3d(nn,i,j,:) = tleaf_shade
                                    tleaf_ave240_tmp_3d(nn,i,j,:)   = tleaf_ave
                                    tmp2mref_tmp_3d(nn,i,j)         = tmp2mref
                                    ubzref_tmp_3d(nn,i,j)           = ubzref

                                    !!TODO:  Restart capability needed to get past leaf temp and PAR if avaialble
                                    !For now, if <= 24 hours then only option is to use current instantaneous values
                                    ppfd_sun24_3d(i,j,:)     = ppfd_sun
                                    ppfd_shade24_3d(i,j,:)   = ppfd_shade
                                    tleaf_sun24_3d(i,j,:)    = tleaf_sun
                                    tleaf_shade24_3d(i,j,:)  = tleaf_shade
                                    tleaf_ave24_3d(i,j,:)    = tleaf_ave
                                    ppfd_sun240_3d(i,j,:)    = ppfd_sun
                                    ppfd_shade240_3d(i,j,:)  = ppfd_shade
                                    tleaf_sun240_3d(i,j,:)   = tleaf_sun
                                    tleaf_shade240_3d(i,j,:) = tleaf_shade
                                    tleaf_ave240_3d(i,j,:)   = tleaf_ave
                                    !for AQ, temp, and wind stress factors in biogenics, still use instantaneous only (caution!)
                                    daily_maxt2m_2d(i,j)        = tmp2mref
                                    daily_mint2m_2d(i,j)        = tmp2mref
                                    daily_maxws10m_2d(i,j)      = ubzref

                                else if (nn .le. 240) then
                                    !Drop nth value into end time slice
                                    ppfd_sun24_tmp_3d(25,i,j,:)     = ppfd_sun
                                    ppfd_shade24_tmp_3d(25,i,j,:)   = ppfd_shade
                                    tleaf_sun24_tmp_3d(25,i,j,:)    = tleaf_sun
                                    tleaf_shade24_tmp_3d(25,i,j,:)  = tleaf_shade
                                    tleaf_ave24_tmp_3d(25,i,j,:)    = tleaf_ave
                                    ppfd_sun240_tmp_3d(nn,i,j,:)    = ppfd_sun
                                    ppfd_shade240_tmp_3d(nn,i,j,:)  = ppfd_shade
                                    tleaf_sun240_tmp_3d(nn,i,j,:)   = tleaf_sun
                                    tleaf_shade240_tmp_3d(nn,i,j,:) = tleaf_shade
                                    tleaf_ave240_tmp_3d(nn,i,j,:)   = tleaf_ave
                                    tmp2mref_tmp_3d(25,i,j)         = tmp2mref
                                    ubzref_tmp_3d(25,i,j)           = ubzref

                                    !Calculate 24 and 240 hour averages
                                    ppfd_sun24_3d(i,j,:)     = sum(ppfd_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_shade24_3d(i,j,:)   = sum(ppfd_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_sun24_3d(i,j,:)    = sum(tleaf_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_shade24_3d(i,j,:)  = sum(tleaf_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_ave24_3d(i,j,:)    = sum(tleaf_ave24_tmp_3d(1:24,i,j,:),1)/24.0
                                    !Use 24 hour data until 240 data become available
                                    ppfd_sun240_3d(i,j,:)    = sum(ppfd_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_shade240_3d(i,j,:)  = sum(ppfd_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_sun240_3d(i,j,:)   = sum(tleaf_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_shade240_3d(i,j,:) = sum(tleaf_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_ave240_3d(i,j,:)   = sum(tleaf_ave24_tmp_3d(1:24,i,j,:),1)/24.0
                                    !Updated rolling 24 hour (hourly, short term) and 240 hour (daily, long-term) averages
                                    ppfd_sun24_3d(i,j,:)     = ( holdfrac * ppfd_sun24_3d(i,j,:) ) + ( hnewfrac * ppfd_sun )
                                    ppfd_shade24_3d(i,j,:)   = ( holdfrac * ppfd_shade24_3d(i,j,:) ) + ( hnewfrac * ppfd_shade )
                                    tleaf_sun24_3d(i,j,:)    = ( holdfrac * tleaf_sun24_3d(i,j,:) ) + ( hnewfrac * tleaf_sun )
                                    tleaf_shade24_3d(i,j,:)  = ( holdfrac * tleaf_shade24_3d(i,j,:) ) + ( hnewfrac * tleaf_shade )
                                    tleaf_ave24_3d(i,j,:)    = ( holdfrac * tleaf_ave24_3d(i,j,:) ) + ( hnewfrac * tleaf_ave )
                                    ppfd_sun240_3d(i,j,:)     = ( holdfrac * ppfd_sun24_3d(i,j,:) ) + ( hnewfrac * ppfd_sun )
                                    ppfd_shade240_3d(i,j,:)   = ( holdfrac * ppfd_shade24_3d(i,j,:) ) + ( hnewfrac * ppfd_shade )
                                    tleaf_sun240_3d(i,j,:)    = ( holdfrac * tleaf_sun24_3d(i,j,:) ) + ( hnewfrac * tleaf_sun )
                                    tleaf_shade240_3d(i,j,:)  = ( holdfrac * tleaf_shade24_3d(i,j,:) ) + ( hnewfrac * tleaf_shade )
                                    tleaf_ave240_3d(i,j,:)    = ( holdfrac * tleaf_ave24_3d(i,j,:) ) + ( hnewfrac * tleaf_ave )
                                    !Take daily max/min for temperature and daily max for wind speed
                                    daily_maxt2m_2d(i,j)        = maxval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_mint2m_2d(i,j)        = minval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_maxws10m_2d(i,j)      = maxval(ubzref_tmp_3d(1:24,i,j))

                                    !Shift time window to make room for new nth value
                                    ppfd_sun24_tmp_3d(1:24,i,j,:)     = ppfd_sun24_tmp_3d(2:25,i,j,:)
                                    ppfd_shade24_tmp_3d(1:24,i,j,:)   = ppfd_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_sun24_tmp_3d(1:24,i,j,:)    = tleaf_sun24_tmp_3d(2:25,i,j,:)
                                    tleaf_shade24_tmp_3d(1:24,i,j,:)  = tleaf_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_ave24_tmp_3d(1:24,i,j,:)    = tleaf_ave24_tmp_3d(2:25,i,j,:)
                                    tmp2mref_tmp_3d(1:24,i,j)         = tmp2mref_tmp_3d(2:25,i,j)
                                    ubzref_tmp_3d(1:24,i,j)           = ubzref_tmp_3d(2:25,i,j)

                                else
                                    !Drop nth value into end time slice
                                    ppfd_sun24_tmp_3d(25,i,j,:)     = ppfd_sun
                                    ppfd_shade24_tmp_3d(25,i,j,:)   = ppfd_shade
                                    tleaf_sun24_tmp_3d(25,i,j,:)    = tleaf_sun
                                    tleaf_shade24_tmp_3d(25,i,j,:)  = tleaf_shade
                                    tleaf_ave24_tmp_3d(25,i,j,:)    = tleaf_ave
                                    ppfd_sun240_tmp_3d(241,i,j,:)    = ppfd_sun
                                    ppfd_shade240_tmp_3d(241,i,j,:)  = ppfd_shade
                                    tleaf_sun240_tmp_3d(241,i,j,:)   = tleaf_sun
                                    tleaf_shade240_tmp_3d(241,i,j,:) = tleaf_shade
                                    tleaf_ave240_tmp_3d(241,i,j,:)   = tleaf_ave
                                    tmp2mref_tmp_3d(25,i,j)         = tmp2mref
                                    ubzref_tmp_3d(25,i,j)           = ubzref

                                    !Calculate 24 and 240 hour averages
                                    ppfd_sun24_3d(i,j,:)     = sum(ppfd_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_shade24_3d(i,j,:)   = sum(ppfd_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_sun24_3d(i,j,:)    = sum(tleaf_sun24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_shade24_3d(i,j,:)  = sum(tleaf_shade24_tmp_3d(1:24,i,j,:),1)/24.0
                                    tleaf_ave24_3d(i,j,:)    = sum(tleaf_ave24_tmp_3d(1:24,i,j,:),1)/24.0
                                    ppfd_sun240_3d(i,j,:)    = sum(ppfd_sun240_tmp_3d(1:240,i,j,:),1)/240.0
                                    ppfd_shade240_3d(i,j,:)  = sum(ppfd_shade240_tmp_3d(1:240,i,j,:),1)/240.0
                                    tleaf_sun240_3d(i,j,:)   = sum(tleaf_sun240_tmp_3d(1:240,i,j,:),1)/240.0
                                    tleaf_shade240_3d(i,j,:) = sum(tleaf_shade240_tmp_3d(1:240,i,j,:),1)/240.0
                                    tleaf_ave240_3d(i,j,:)   = sum(tleaf_ave240_tmp_3d(1:240,i,j,:),1)/240.0
                                    !Update for current time value using efolding (holdfrac and hnewfrac)
                                    ppfd_sun24_3d(i,j,:)     = ( holdfrac * ppfd_sun24_3d(i,j,:) ) + ( hnewfrac * ppfd_sun )
                                    ppfd_shade24_3d(i,j,:)   = ( holdfrac * ppfd_shade24_3d(i,j,:) ) + ( hnewfrac * ppfd_shade )
                                    tleaf_sun24_3d(i,j,:)    = ( holdfrac * tleaf_sun24_3d(i,j,:) ) + ( hnewfrac * tleaf_sun )
                                    tleaf_shade24_3d(i,j,:)  = ( holdfrac * tleaf_shade24_3d(i,j,:) ) + ( hnewfrac * tleaf_shade )
                                    tleaf_ave24_3d(i,j,:)    = ( holdfrac * tleaf_ave24_3d(i,j,:) ) + ( hnewfrac * tleaf_ave )
                                    ppfd_sun240_3d(i,j,:)    = ( doldfrac * ppfd_sun240_3d(i,j,:) ) + ( dnewfrac * ppfd_sun )
                                    ppfd_shade240_3d(i,j,:)  = ( doldfrac * ppfd_shade240_3d(i,j,:) ) + ( dnewfrac * ppfd_shade )
                                    tleaf_sun240_3d(i,j,:)   = ( doldfrac * tleaf_sun240_3d(i,j,:) ) + ( dnewfrac * tleaf_sun )
                                    tleaf_shade240_3d(i,j,:) = ( doldfrac * tleaf_shade240_3d(i,j,:) ) + ( dnewfrac * tleaf_shade )
                                    tleaf_ave240_3d(i,j,:)   = ( doldfrac * tleaf_ave240_3d(i,j,:) ) + ( dnewfrac * tleaf_ave )
                                    !Take daily max/min for temperature and daily max for wind speed
                                    daily_maxt2m_2d(i,j)        = maxval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_mint2m_2d(i,j)        = minval(tmp2mref_tmp_3d(1:24,i,j))
                                    daily_maxws10m_2d(i,j)      = maxval(ubzref_tmp_3d(1:24,i,j))

                                    !Shift time window to make room for new nth value
                                    ppfd_sun24_tmp_3d(1:24,i,j,:)     = ppfd_sun24_tmp_3d(2:25,i,j,:)
                                    ppfd_shade24_tmp_3d(1:24,i,j,:)   = ppfd_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_sun24_tmp_3d(1:24,i,j,:)    = tleaf_sun24_tmp_3d(2:25,i,j,:)
                                    tleaf_shade24_tmp_3d(1:24,i,j,:)  = tleaf_shade24_tmp_3d(2:25,i,j,:)
                                    tleaf_ave24_tmp_3d(1:24,i,j,:)    = tleaf_ave24_tmp_3d(2:25,i,j,:)
                                    ppfd_sun240_tmp_3d(1:240,i,j,:)    = ppfd_sun240_tmp_3d(2:241,i,j,:)
                                    ppfd_shade240_tmp_3d(1:240,i,j,:)  = ppfd_shade240_tmp_3d(2:241,i,j,:)
                                    tleaf_sun240_tmp_3d(1:240,i,j,:)   = tleaf_sun240_tmp_3d(2:241,i,j,:)
                                    tleaf_shade240_tmp_3d(1:240,i,j,:) = tleaf_shade240_tmp_3d(2:241,i,j,:)
                                    tleaf_ave240_tmp_3d(1:240,i,j,:)   = tleaf_ave240_tmp_3d(2:241,i,j,:)
                                    tmp2mref_tmp_3d(1:24,i,j)         = tmp2mref_tmp_3d(2:25,i,j)
                                    ubzref_tmp_3d(1:24,i,j)           = ubzref_tmp_3d(2:25,i,j)
                                end if
                            else
                                write(*,*) 'wrong HIST_OPT namelist option = ', hist_opt, 'only option = 0 or 1'
                                call exit(2)
                            end if
! ... user option to calculate in-canopy biogenic emissions
                            if (ifcanbio) then
                                if (cluref .gt. 0.0_rk) then
                                    !ISOP
                                    if (biospec_opt == 0 .or. biospec_opt == 1) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 1, emi_isop_3d(i,j,:))
                                    end if
                                    !MYRC
                                    if (biospec_opt == 0 .or. biospec_opt == 2) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 2, emi_myrc_3d(i,j,:))
                                    end if
                                    !SABI
                                    if (biospec_opt == 0 .or. biospec_opt == 3) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 3, emi_sabi_3d(i,j,:))
                                    end if
                                    !LIMO
                                    if (biospec_opt == 0 .or. biospec_opt == 4) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 4, emi_limo_3d(i,j,:))
                                    end if
                                    !CARE
                                    if (biospec_opt == 0 .or. biospec_opt == 5) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 5, emi_care_3d(i,j,:))
                                    end if
                                    !OCIM
                                    if (biospec_opt == 0 .or. biospec_opt == 6) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 6, emi_ocim_3d(i,j,:))
                                    end if
                                    !BPIN
                                    if (biospec_opt == 0 .or. biospec_opt == 7) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 7, emi_bpin_3d(i,j,:))
                                    end if
                                    !APIN
                                    if (biospec_opt == 0 .or. biospec_opt == 8) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 8, emi_apin_3d(i,j,:))
                                    end if
                                    !MONO
                                    if (biospec_opt == 0 .or. biospec_opt == 9) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 9, emi_mono_3d(i,j,:))
                                    end if
                                    !FARN
                                    if (biospec_opt == 0 .or. biospec_opt == 10) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 10, emi_farn_3d(i,j,:))
                                    end if
                                    !CARY
                                    if (biospec_opt == 0 .or. biospec_opt == 11) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 11, emi_cary_3d(i,j,:))
                                    end if
                                    !SESQ
                                    if (biospec_opt == 0 .or. biospec_opt == 12) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 12, emi_sesq_3d(i,j,:))
                                    end if
                                    !MBOL
                                    if (biospec_opt == 0 .or. biospec_opt == 13) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 13, emi_mbol_3d(i,j,:))
                                    end if
                                    !METH
                                    if (biospec_opt == 0 .or. biospec_opt == 14) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 14, emi_meth_3d(i,j,:))
                                    end if
                                    !ACET
                                    if (biospec_opt == 0 .or. biospec_opt == 15) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 15, emi_acet_3d(i,j,:))
                                    end if
                                    !CO
                                    if (biospec_opt == 0 .or. biospec_opt == 16) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 16, emi_co_3d(i,j,:))
                                    end if
                                    !BIDI VOC
                                    if (biospec_opt == 0 .or. biospec_opt == 17) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 17, emi_bvoc_3d(i,j,:))
                                    end if
                                    !Stress VOC
                                    if (biospec_opt == 0 .or. biospec_opt == 18) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 18, emi_svoc_3d(i,j,:))
                                    end if
                                    !Other VOC
                                    if (biospec_opt == 0 .or. biospec_opt == 19) then
                                        call canopy_bio(zk, fafraczint, hcmref, &
                                            lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                            ppfd_sun24_3d(i,j,:), ppfd_shade24_3d(i,j,:), &
                                            tleaf_ave24_3d(i,j,:), ppfd_sun240_3d(i,j,:), ppfd_shade240_3d(i,j,:), &
                                            tleaf_ave240_3d(i,j,:), tka_3d(i,j,:), dswrfref, tmp2mref, &
                                            lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                            leafage_opt, pastlai, currentlai, tsteplai,  &
                                            loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                            soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                            soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                            ht_opt, lt_opt, hw_opt, daily_maxt2m_2d(i,j), daily_mint2m_2d(i,j), &
                                            daily_maxws10m_2d(i,j), &
                                            modlays, 19, emi_ovoc_3d(i,j,:))
                                    end if
                                else
                                    if (biospec_opt == 0 .or. biospec_opt == 1) then
                                        emi_isop_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 2) then
                                        emi_myrc_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 3) then
                                        emi_sabi_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 4) then
                                        emi_limo_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 5) then
                                        emi_care_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 6) then
                                        emi_ocim_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 7) then
                                        emi_bpin_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 8) then
                                        emi_apin_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 9) then
                                        emi_mono_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 10) then
                                        emi_farn_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 11) then
                                        emi_cary_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 12) then
                                        emi_sesq_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 13) then
                                        emi_mbol_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 14) then
                                        emi_meth_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 15) then
                                        emi_acet_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 16) then
                                        emi_co_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 17) then
                                        emi_bvoc_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 18) then
                                        emi_svoc_3d(i,j,:) = 0.0_rk
                                    end if
                                    if (biospec_opt == 0 .or. biospec_opt == 19) then
                                        emi_ovoc_3d(i,j,:) = 0.0_rk
                                    end if
                                end if
                            end if

! ... user option to calculate in-canopy dry deposition velocity
                            if (ifcanddepgas ) then
                                if (ifcanwind) then   !ubar needed for rbl
                                    if (chemmechgas_opt .eq. 0)  then   !RACM2
                                        if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                            !At this point must have canopy height set by thresholds, so Ra calculation is valid
                                            !Calculate Bulk Richardson Number for stability dependent aerodynamic resistance
                                            rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                                d_h_2d(i,j)*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk)
                                            ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, d_h_2d(i,j)*hcmref*100.0_rk, &
                                                hcmref*100.0_rk, rib)
                                            ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum,and assume constant Ra
                                            !downward at all levels through canopy to ground
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra, 1, ddep_no_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra, 1, ddep_no_3d(i,j,1))   ! [cm/s]
                                                else !snow or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 1,ddep_no_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  2, ddep_no2_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra, 2, ddep_no2_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 2,ddep_no2_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 3, ddep_o3_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref,  ra, 3, ddep_o3_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 3,ddep_o3_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 4, ddep_hono_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref,  ra, 4, ddep_hono_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 4,ddep_hono_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 5, ddep_hno4_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref,  ra, 5, ddep_hno4_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 5,ddep_hno4_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  6, ddep_hno3_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  6, ddep_hno3_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 6,ddep_hno3_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 7, ddep_n2o5_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  7, ddep_n2o5_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 7,ddep_n2o5_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  8, ddep_co_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref,  ra, 8, ddep_co_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 8,ddep_co_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 9, ddep_h2o2_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref,  ra, 9, ddep_h2o2_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 9,ddep_h2o2_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 10, ddep_ch4_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref,  ra, 10, ddep_ch4_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 10,ddep_ch4_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 11, ddep_mo2_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  11, ddep_mo2_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 11,ddep_mo2_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  12, ddep_op1_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,2), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  12, ddep_op1_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 12,ddep_op1_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  13, ddep_moh_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  13, ddep_moh_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 13,ddep_moh_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref,  ra, 14, ddep_no3_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  14, ddep_no3_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 14,ddep_no3_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  15, ddep_o3p_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  15, ddep_o3p_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 15,ddep_o3p_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  16, ddep_o1d_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  16, ddep_o1d_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 16,ddep_o1d_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  17, ddep_ho_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  17, ddep_ho_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 17,ddep_ho_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  18, ddep_ho2_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  18, ddep_ho2_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 18,ddep_ho2_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  19, ddep_ora1_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref,  ra, 19, ddep_ora1_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 19,ddep_ora1_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  20, ddep_hac_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  20, ddep_hac_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 20,ddep_hac_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  21, ddep_paa_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  21, ddep_paa_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 21,ddep_paa_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  22, ddep_dhmob_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  22, ddep_dhmob_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 22,ddep_dhmob_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  23, ddep_hpald_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  23, ddep_hpald_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 23,ddep_hpald_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  24, ddep_ishp_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  24, ddep_ishp_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 24,ddep_ishp_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  25, ddep_iepox_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  25, ddep_iepox_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 25,ddep_iepox_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  26, ddep_propnn_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  26, ddep_propnn_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 26,ddep_propnn_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  27, ddep_isopnb_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  27, ddep_isopnb_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 27,ddep_isopnb_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  28, ddep_isopnd_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  28, ddep_isopnd_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 28,ddep_isopnd_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  29, ddep_macrn_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  29, ddep_macrn_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 29,ddep_macrn_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  30, ddep_mvkn_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  30, ddep_mvkn_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 30,ddep_mvkn_3d(i,j,1))
                                                end if
                                            endif
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                                call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                    zk, hcmref, tka_3d(i,j,:), pressa_3d(i,j,:), &
                                                    relhuma_3d(i,j,:), fsun, ppfd_sun, ppfd_shade, canwind_3d(i,j,:),  &
                                                    dswrfref, ra,  31, ddep_isnp_3d(i,j,:))   ! [cm/s]
                                            endif
                                            !soil gas dry depostion at level 1, i.e., z=0
                                            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                                if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                    call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                        soilt1ref, pressa_3d(i,j,1), canwind_3d(i,j,2),  &
                                                        soilcat_opt,sotypref,soild1, soilw1ref, ra,  31, ddep_isnp_3d(i,j,1))   ! [cm/s]
                                                else !snow covered or ice covered
                                                    call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                        canwind_3d(i,j,2), ra, 31,ddep_isnp_3d(i,j,1))
                                                end if
                                            endif
                                        else
                                            write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                            write(*,*)  ' in namelist...exiting'
                                            write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                            call exit(2)
                                        end if
                                    else
                                        write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                                        write(*,*)  'Set chemmechgas_opt = 0 (RACM2) for now'
                                        call exit(2)
                                    end if
                                else
                                    write(*,*)  'Wrong IfCanWind choice of ', ifcanwind, ' in namelist...exiting'
                                    write(*,*)  'Set IfCanwind to True to use IfCanDDepGas'
                                    call exit(2)
                                end if
                            end if

                            ! --- Sub-canopy aerosol dry deposition through vegetative canopies (Katul et al. 2010) ---
                            if (ifcanaeroddep) then
                                if (ifcanwind) then !ubar needed for rbl
                                    rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                                        d_h_2d(i,j)*hcmref, (hcmref+hgtref))
                                    ra=rav(ubzref, (hcmref+hgtref), d_h_2d(i,j)*hcmref, &
                                        hcmref, rib)
                                    ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                                    call canopy_aero_ddep_subveg(modlays, zk, hcmref, lad_3d(i,j,:), canwind_3d(i,j,:), &
                                        aeroddep_diam, aeroddep_rho, tka_3d(i,j,:), pressa_3d(i,j,:), aeroddep_opt, &
                                        ra, modres, ustref, vtyperef, vdep_aero_3d(i,j,:)) ! [cm/s]
                                else
                                    write(*,*)  'Wrong IfCanWind choice of ', ifcanwind, ' in namelist...exiting'
                                    write(*,*)  'Set IfCanwind to True to use IfCanAeroDDep'
                                    call exit(2)
                                end if
                                !  --- Sub-canopy aerosol dry deposition to soil (i.e., first model layer) beneath canopy (Pleim et al. 2022) ---
                                rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                                    d_h_2d(i,j)*hcmref, (hcmref+hgtref))
                                ra=rav(ubzref, (hcmref+hgtref), d_h_2d(i,j)*hcmref, &
                                    hcmref, rib)
                                ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                                call canopy_aero_ddep_pleim2022(ustref, ra, canwind_3d(i,j,2), &
                                    aeroddep_diam, aeroddep_rho, tka_3d(i,j,1), pressa_3d(i,j,1), 2, vdep_aero_3d(i,j,1))
                            end if

                        else
                            if (biospec_opt == 0 .or. biospec_opt == 1) then
                                emi_isop_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 2) then
                                emi_myrc_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 3) then
                                emi_sabi_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 4) then
                                emi_limo_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 5) then
                                emi_care_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 6) then
                                emi_ocim_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 7) then
                                emi_bpin_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 8) then
                                emi_apin_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 9) then
                                emi_mono_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 10) then
                                emi_farn_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 11) then
                                emi_cary_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 12) then
                                emi_sesq_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 13) then
                                emi_mbol_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 14) then
                                emi_meth_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 15) then
                                emi_acet_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 16) then
                                emi_co_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 17) then
                                emi_bvoc_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 18) then
                                emi_svoc_3d(i,j,:) = 0.0_rk
                            end if
                            if (biospec_opt == 0 .or. biospec_opt == 19) then
                                emi_ovoc_3d(i,j,:) = 0.0_rk
                            end if
                        end if !Contiguous Canopy

                    else if (vtyperef .eq. 15 .or. vtyperef .eq. 16 .or. vtyperef .eq. 20) then !Barren/Sparsely Vegetated  or
                        !Barren Tundra or Snow/Ice
! ... user option to calculate dry deposition velocity...for land use outside of vegetated canopies
                        if (ifcanddepgas ) then
                            if (chemmechgas_opt .eq. 0)  then   !RACM2
                                if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (hcmref .le. 0.0) then
                                        hcmref = 0.01!m   (set to just above ground @ 1 cm; leads to z0 ~ 0.001 m
                                    endif                ! needed for Ra calculation)
                                    rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                        0.75_rk*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk) !assume d~3/4*Hc
                                    ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, 0.75_rk*hcmref*100.0_rk, &
                                        hcmref*100.0_rk, rib)
                                    ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  1, ddep_no_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 1,ddep_no_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  2, ddep_no2_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 2,ddep_no2_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  3, ddep_o3_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 3,ddep_o3_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  4, ddep_hono_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 4,ddep_hono_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  5, ddep_hno4_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 5,ddep_hno4_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  6, ddep_hno3_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 6,ddep_hno3_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  7, ddep_n2o5_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 7,ddep_n2o5_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  8, ddep_co_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 8,ddep_co_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  9, ddep_h2o2_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 9,ddep_h2o2_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  10, ddep_ch4_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 10,ddep_ch4_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  11, ddep_mo2_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 11,ddep_mo2_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  12, ddep_op1_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 12,ddep_op1_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  13, ddep_moh_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 13,ddep_moh_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  14, ddep_no3_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 14,ddep_no3_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  15, ddep_o3p_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 15,ddep_o3p_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  16, ddep_o1d_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 16,ddep_o1d_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  17, ddep_ho_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 17,ddep_ho_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  18, ddep_ho2_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 18,ddep_ho2_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  19, ddep_ora1_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 19,ddep_ora1_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  20, ddep_hac_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 20,ddep_hac_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  21, ddep_paa_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 21,ddep_paa_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  22, ddep_dhmob_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 22,ddep_dhmob_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  23, ddep_hpald_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 23,ddep_hpald_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  24, ddep_ishp_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 24,ddep_ishp_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  25, ddep_iepox_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 25,ddep_iepox_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  26, ddep_propnn_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 26,ddep_propnn_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  27, ddep_isopnb_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 27,ddep_isopnb_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  28, ddep_isopnd_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 28,ddep_isopnd_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  29, ddep_macrn_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 29,ddep_macrn_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  30, ddep_mvkn_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 30,ddep_mvkn_3d(i,j,1))
                                        endif
                                    endif
                                    !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                soilt1ref, pressa_3d(i,j,1), ubzref,  &
                                                soilcat_opt,sotypref,soild1, soilw1ref, ra,  31, ddep_isnp_3d(i,j,1))   ! [cm/s]
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 31,ddep_isnp_3d(i,j,1))
                                        endif
                                    endif
                                else
                                    write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                    write(*,*)  ' in namelist...exiting'
                                    write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                    call exit(2)
                                end if
                            else
                                write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                                write(*,*)  'Set chemmechgas_opt = 0 (RACM2) for now'
                                call exit(2)
                            end if
                        end if
                        ! --- Sub-canopy aerosol dry deposition over barren/soil surfaces (Pleim et al. 2022) ---
                        if (ifcanaeroddep) then
                            rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                                d_h_2d(i,j)*hcmref, (hcmref+hgtref))
                            ra=rav(ubzref, (hcmref+hgtref), d_h_2d(i,j)*hcmref, &
                                hcmref, rib)
                            ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                            call canopy_aero_ddep_pleim2022(ustref, ra, ubzref, &
                                aeroddep_diam, aeroddep_rho, tmp2mref, pressfcref, 2, vdep_aero_3d(i,j,1))
                        end if
                    else if (vtyperef .eq. 13) then !Urban and Built Up
! ... user option to calculate dry deposition velocity...for land use outside of vegetated canopies
                        if (ifcanddepgas ) then
                            if (chemmechgas_opt .eq. 0)  then   !RACM2
                                if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                    !urban gas dry depostion at level 1, i.e., z=0
                                    if (hcmref .le. 0.0) then
                                        hcmref = 0.01!m   (set to just above ground @ 1 cm; leads to z0 ~ 0.001 m
                                    endif                ! needed for Ra calculation)
                                    rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                        0.75_rk*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk) !assume d~3/4*Hc
                                    ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, 0.75_rk*hcmref*100.0_rk, &
                                        hcmref*100.0_rk, rib)
                                    ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,1,ddep_no_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 1,ddep_no_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,2,ddep_no2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 2,ddep_no2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,3,ddep_o3_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 3,ddep_o3_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,4,ddep_hono_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 4,ddep_hono_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,5,ddep_hno4_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 5,ddep_hno4_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,6,ddep_hno3_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 6,ddep_hno3_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,7,ddep_n2o5_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 7,ddep_n2o5_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,8,ddep_co_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 8,ddep_co_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,9,ddep_h2o2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 9,ddep_h2o2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,10,ddep_ch4_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 10,ddep_ch4_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,11,ddep_mo2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 11,ddep_mo2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,12,ddep_op1_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 12,ddep_op1_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,13,ddep_moh_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 13,ddep_moh_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,14,ddep_no3_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 14,ddep_no3_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,15,ddep_o3p_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 15,ddep_o3p_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,16,ddep_o1d_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 16,ddep_o1d_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,17,ddep_ho_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 17,ddep_ho_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,18,ddep_ho2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 18,ddep_ho2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,19,ddep_ora1_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 19,ddep_ora1_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,20,ddep_hac_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 20,ddep_hac_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,21,ddep_paa_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 21,ddep_paa_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,22,ddep_dhmob_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 22,ddep_dhmob_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,23,ddep_hpald_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 23,ddep_hpald_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,24,ddep_ishp_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 24,ddep_ishp_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,25,ddep_iepox_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 25,ddep_iepox_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,26,ddep_propnn_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 26,ddep_propnn_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,27,ddep_isopnb_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 27,ddep_isopnb_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,28,ddep_isopnd_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 28,ddep_isopnd_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,29,ddep_macrn_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 29,ddep_macrn_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,30,ddep_mvkn_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 30,ddep_mvkn_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,tmp2mref,gamma_set, ra,31,ddep_isnp_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 31,ddep_isnp_3d(i,j,1))
                                        endif
                                    endif
                                else
                                    write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                    write(*,*)  ' in namelist...exiting'
                                    write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                    call exit(2)
                                end if
                            else
                                write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                                write(*,*)  'Set chemmechgas_opt = 0 (RACM2) for now'
                                call exit(2)
                            end if
                        end if
                        ! --- Sub-canopy aerosol dry deposition over urban surfaces (Pleim et al. 2022) ---
                        if (ifcanaeroddep) then
                            rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                                d_h_2d(i,j)*hcmref, (hcmref+hgtref))
                            ra=rav(ubzref, (hcmref+hgtref), d_h_2d(i,j)*hcmref, &
                                hcmref, rib)
                            ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                            call canopy_aero_ddep_pleim2022(ustref, ra, ubzref, &
                                aeroddep_diam, aeroddep_rho, tmp2mref, pressfcref, 1, vdep_aero_3d(i,j,1))
                        end if
                    else if (vtyperef .eq. 0) then !Water from FV3 (usually vtype = 17 for water)
! ... user option to calculate dry deposition velocity...for land use outside of vegetated canopies
                        if (ifcanddepgas ) then
                            if (chemmechgas_opt .eq. 0)  then   !RACM2
                                if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                    !water gas dry depostion at level 1, i.e., z=0
                                    if (hcmref .le. 0.0) then
                                        hcmref = 0.01!m   (set to just above ground @ 1 cm; leads to z0 ~ 0.001 m
                                    endif                ! needed for Ra calculation)
                                    rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                        0.75_rk*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk) !assume d~3/4*Hc
                                    ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, 0.75_rk*hcmref*100.0_rk, &
                                        hcmref*100.0_rk, rib)
                                    ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,1,ddep_no_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 1,ddep_no_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,2,ddep_no2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 2,ddep_no2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,3,ddep_o3_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 3,ddep_o3_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,4,ddep_hono_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 4,ddep_hono_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,5,ddep_hno4_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 5,ddep_hno4_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,6,ddep_hno3_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 6,ddep_hno3_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,7,ddep_n2o5_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 7,ddep_n2o5_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,8,ddep_co_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 8,ddep_co_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,9,ddep_h2o2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 9,ddep_h2o2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,10,ddep_ch4_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 10,ddep_ch4_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,11,ddep_mo2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 11,ddep_mo2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,12,ddep_op1_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 12,ddep_op1_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,13,ddep_moh_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 13,ddep_moh_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,14,ddep_no3_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 14,ddep_no3_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,15,ddep_o3p_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 15,ddep_o3p_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,16,ddep_o1d_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 16,ddep_o1d_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,17,ddep_ho_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 17,ddep_ho_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,18,ddep_ho2_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref,  ra,18,ddep_ho2_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,19,ddep_ora1_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 19,ddep_ora1_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,20,ddep_hac_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 20,ddep_hac_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,21,ddep_paa_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 21,ddep_paa_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,22,ddep_dhmob_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 22,ddep_dhmob_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,23,ddep_hpald_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 23,ddep_hpald_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,24,ddep_ishp_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 24,ddep_ishp_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,25,ddep_iepox_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 25,ddep_iepox_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,26,ddep_propnn_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 26,ddep_propnn_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,27,ddep_isopnb_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 27,ddep_isopnb_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,28,ddep_isopnd_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 28,ddep_isopnd_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,29,ddep_macrn_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 29,ddep_macrn_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,30,ddep_mvkn_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 30,ddep_mvkn_3d(i,j,1))
                                        endif
                                    endif
                                    if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                        if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                            call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                                tmp2mref,spfh2mref,ustref, ra,31,ddep_isnp_3d(i,j,1))
                                        else  !snow covered or ice covered
                                            call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                ubzref, ra, 31,ddep_isnp_3d(i,j,1))
                                        endif
                                    endif
                                else
                                    write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                    write(*,*)  ' in namelist...exiting'
                                    write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                    call exit(2)
                                end if
                            else
                                write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                                write(*,*)  'Set chemmechgas_opt = 0 (RACM2) for now'
                                call exit(2)
                            end if
                        end if
                        ! --- Sub-canopy aerosol dry deposition over default surfaces including water (Pleim et al. 2022) ---
                        if (ifcanaeroddep) then
                            rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                                d_h_2d(i,j)*hcmref, (hcmref+hgtref))
                            ra=rav(ubzref, (hcmref+hgtref), d_h_2d(i,j)*hcmref, &
                                hcmref, rib)
                            ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                            call canopy_aero_ddep_pleim2022(ustref, ra, ubzref, &
                                aeroddep_diam, aeroddep_rho, tmp2mref, pressfcref, 3, vdep_aero_3d(i,j,1))
                        end if
                    else
                        write(*,*)  'Warning VIIRS/MODIS VTYPE ', vtyperef, ' is not supported...continue'
                    end if   !Vegetation types
                else
                    write(*,*)  'Wrong LU_OPT choice of ', lu_opt, ' in namelist...exiting'
                    call exit(2)

                end if       !Landuse Options
            end do !lat
        end do   !lon

!----------------------------------------------------------->

    else if (infmt_opt .eq. 1) then !Main input format is 1D/2D text and output will be 1D/2D text

        if (ifcanwind .or. ifcanwaf) then !only calculate if canopy wind or WAF option
            call canopy_calcdx(dx_opt, dx_set, nlat, nlon, lat1d, &
                lon1d, dx)
        end if

        if (href_opt .eq. 0 ) then !setting entire array = href_set value from user NL
            variables%href = href_set
        else if (href_opt .eq. 1 ) then !from file array
            variables%href =  variables%href
        else
            write(*,*)  'Wrong HREF_OPT choice of ', href_opt, ' in namelist...exiting'
            call exit(2)
        end if

! ... Main loop through model grid cells
        do loc=1, nlat*nlon
            hcmref         = variables(loc)%ch
            uref           = variables(loc)%ugrd10m
            vref           = variables(loc)%vgrd10m
            cluref         = variables(loc)%clu
            lairef         = variables(loc)%lai
            vtyperef       = variables(loc)%vtype
            canfracref     = variables(loc)%canfrac
            ustref         = variables(loc)%fricv
            cszref         = variables(loc)%csz
            z0ref          = variables(loc)%sfcr
            molref         = variables(loc)%mol
            frpref         = variables(loc)%frp
            hgtref         = variables(loc)%href
            sotypref       = variables(loc)%sotyp
            pressfcref     = variables(loc)%pressfc
            dswrfref       = variables(loc)%dswrf
            shtflref       = variables(loc)%shtfl
            tmpsfcref      = variables(loc)%tmpsfc
            tmp2mref       = variables(loc)%tmp2m
            spfh2mref      = variables(loc)%spfh2m
            hpblref        = variables(loc)%hpbl
            prate_averef   = variables(loc)%prate_ave
            soilw1ref      = variables(loc)%soilw1
            soilw2ref      = variables(loc)%soilw2
            soilw3ref      = variables(loc)%soilw3
            soilw4ref      = variables(loc)%soilw4
            wiltref        = variables(loc)%wilt
            ozone_w126ref  = variables(loc)%ozone_w126
            soilt1ref      = variables(loc)%soilt1
            soilt2ref      = variables(loc)%soilt2
            soilt3ref      = variables(loc)%soilt3
            soilt4ref      = variables(loc)%soilt4
            tmp_hyblev1ref = variables(loc)%tmp_hyblev1
            snowc_averef   = variables(loc)%snowc_ave
            icec_averef    = variables(loc)%icec
            vbdsf_averef   = variables(loc)%vbdsf_ave
            vddsf_averef   = variables(loc)%vddsf_ave

            if (var3d_opt .eq. 1) then !allocated so set
                pavd_arr     = (/variables_can(loc)%pavd01, &
                    variables_can(loc)%pavd02, &
                    variables_can(loc)%pavd03, &
                    variables_can(loc)%pavd04, &
                    variables_can(loc)%pavd05, &
                    variables_can(loc)%pavd06, &
                    variables_can(loc)%pavd07, &
                    variables_can(loc)%pavd08, &
                    variables_can(loc)%pavd09, &
                    variables_can(loc)%pavd10, &
                    variables_can(loc)%pavd11, &
                    variables_can(loc)%pavd12, &
                    variables_can(loc)%pavd13, &
                    variables_can(loc)%pavd14/)
                lev_arr      = (/variables_can(loc)%lev01, &
                    variables_can(loc)%lev02, &
                    variables_can(loc)%lev03, &
                    variables_can(loc)%lev04, &
                    variables_can(loc)%lev05, &
                    variables_can(loc)%lev06, &
                    variables_can(loc)%lev07, &
                    variables_can(loc)%lev08, &
                    variables_can(loc)%lev09, &
                    variables_can(loc)%lev10, &
                    variables_can(loc)%lev11, &
                    variables_can(loc)%lev12, &
                    variables_can(loc)%lev13, &
                    variables_can(loc)%lev14/)
            end if

! ... calculate wind speed from u and v
            ubzref   = sqrt((uref**2.0) + (vref**2.0))

! ... get scaled canopy model profile and sub-canopy layers
            if (hcmref > 0) then
                zhc = zk / hcmref
            else
                zhc = 0
            end if
            cansublays  = min(floor(hcmref/modres),modlays)
            if (cansublays .lt. 1) then !case where model resolution >= canopy height
                cansublays=1            !only one layer allowed
            end if

! ... calculate initial canopy temp/pressure/humidity/density profiles from first guess approximations (i.e., no leaf energy balance)

            do k=1, modlays
                !For temperature we use 1st model layer and 2-m temperatures to estimate lapse rate, not tmpsfc
                !which is too extreme in some global locations.
                tka(loc,k)=calctemp(zk(k)*100.0_rk, (hyblev1*100.0_rk) - 200.0_rk, &
                    tmp_hyblev1ref-273.15_rk, tmp2mref-273.15_rk) ! temp    [K]
                pressa(loc,k)=calcpressure(zk(k)*100.0_rk, 200.0_rk, pressfcref*0.01_rk, &
                    tmp2mref, tmpsfcref)                                                               ! press   [mb]
                relhuma(loc,k)=calcrelhum(tka(loc,k),pressa(loc,k),spfh2mref*1000.0_rk)          ! relhum  [%]
                spechuma(loc,k)=calcspechum(relhuma(loc,k),tka(loc,k),pressa(loc,k))              ! spechum [g/kg]
            end do

! ... check for valid model vegetation types
            if (lu_opt .eq. 0 .or. lu_opt .eq. 1 ) then !VIIRS or MODIS
                if (vtyperef .gt. 0 .and. vtyperef .le. 12 &
                    .or. vtyperef .eq. 14 .or. vtyperef .eq. 18 &
                    .or. vtyperef .eq. 19) then


! ... check for can_opt from user namelist
                    if (vtyperef .ge. 1 .and. vtyperef .le. 5 &
                        .or. vtyperef .eq. 18 .or. vtyperef .eq. 19) then !!VIIRS/MODIS forest/tundra canopies
                        if (can_opt .eq. 0) then !use inputs for canopy heights
                            hcmref = hcmref
                            canfracref = canfracref
                            lairef = lairef
                        else if (can_opt .eq. 1) then !user set constant canopy heights for forests/tundras
                            !only used to overide valid vtypes with hcm<=0, canfrac<=0, or lai<=0.
                            if (hcmref .le. 0.0) then
                                hcmref = can_chset
                            else
                                hcmref = hcmref
                            endif
                            if (canfracref .le. 0.0) then
                                canfracref = can_cfset
                            else
                                canfracref = canfracref
                            endif
                            if (lairef .le. 0.0) then
                                lairef = can_laiset
                            else
                                lairef = lairef
                            endif
                            !recalculate
                            zhc         = zk/hcmref
                            cansublays  = min(floor(hcmref/modres),modlays)
                            if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                cansublays=1            !only one layer allowed
                            end if
                        else
                            write(*,*)  'Wrong CAN_OPT choice of ', can_opt, &
                                ' in namelist...exiting'
                            call exit(2)
                        end if
                    end if

! ... check for ssg_opt from user namelist
                    if (vtyperef .ge. 6 .and. vtyperef .le. 11) then !VIIRS/MODIS shrubs/savannas/grasses (SSG) type
                        !Includes wetlands!
                        if (ssg_opt .eq. 0) then !use GEDI inputs for SSG heights (possibly not measured...)
                            hcmref = hcmref
                            canfracref = canfracref
                            lairef = lairef
                        else if (ssg_opt .eq. 1) then !user set constant shrubs/savannas/grasslands height
                            !only used to overide valid vtypes with hcm<=0, canfrac<=0, or lai<=0.
                            if (hcmref .le. 0.0) then
                                hcmref = ssg_chset
                            else
                                hcmref = hcmref
                            endif
                            if (canfracref .le. 0.0) then
                                canfracref = ssg_cfset
                            else
                                canfracref = canfracref
                            endif
                            if (lairef .le. 0.0) then
                                lairef = ssg_laiset
                            else
                                lairef = lairef
                            endif
                            !recalculate
                            zhc         = zk/hcmref
                            cansublays  = min(floor(hcmref/modres),modlays)
                            if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                cansublays=1            !only one layer allowed
                            end if
                        else
                            write(*,*)  'Wrong SSG_OPT choice of ', ssg_opt, &
                                ' in namelist...exiting'
                            call exit(2)
                        end if
                    end if

! ... check for crop_opt from user namelist
                    if (vtyperef .eq. 12 .or. vtyperef .eq. 14) then !VIIRS/MODIS crop type
                        if (crop_opt .eq. 0) then !use GEDI inputs for crop height
                            hcmref = hcmref
                            canfracref = canfracref
                            lairef = lairef
                        else if (crop_opt .eq. 1) then !user set constant crop height
                            !only used to overide valid vtypes with hcm<=0, canfrac<=0, or lai<=0.
                            if (hcmref .le. 0.0) then
                                hcmref = crop_chset
                            else
                                hcmref = hcmref
                            endif
                            if (canfracref .le. 0.0) then
                                canfracref = crop_cfset
                            else
                                canfracref = canfracref
                            endif
                            if (lairef .le. 0.0) then
                                lairef = crop_laiset
                            else
                                lairef = lairef
                            endif
                            !recalculate
                            zhc         = zk/hcmref
                            cansublays  = min(floor(hcmref/modres),modlays)
                            if (cansublays .lt. 1) then !case where model resolution >= canopy height
                                cansublays=1            !only one layer allowed
                            end if
                        else
                            write(*,*)  'Wrong CROP_OPT choice of ', crop_opt, &
                                ' in namelist...exiting'
                            call exit(2)
                        end if
                    end if

! ... check for contiguous canopy conditions at each model grid cell
                    if (hcmref .gt. ch_thresh .and. canfracref .gt. cf_thresh &
                        .and. lairef .gt. lai_thresh) then

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

                        call canopy_parm(vtyperef, hcmref, canfracref, lairef, &
                            pai_opt, pai_set, lu_opt, firetype, cdrag, &
                            pai, zcanmax, sigmau, sigma1)

! ... Choose between prescribed canopy/foliate shape profile or observed GEDI PAVD profile

                        if (pavd_opt .eq. 0) then
! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
                            call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                fafraczint)
                        else
                            if (var3d_opt .ne. 1) then
                                write(*,*) 'wrong VAR3D_OPT namelist option = ', var3d_opt, &
                                    'change to 1 for supporting PAVD text read'
                                call exit(2)
                            end if
! ... derive canopy/foliage distribution shape profile from interpolated GEDI PAVD profile - bottom up total in-canopy and fraction at z
                            if (variables(loc)%lat .gt. (-1.0_rk*pavd_set) .and. &
                                variables(loc)%lat .lt. pavd_set) then !use GEDI PAVD
                                call canopy_pavd2fafrac(zcanmax, sigmau, sigma1, hcmref, zhc, &
                                    pavd_arr, lev_arr, fafraczint)
                                !check if there is observed canopy height but no PAVD profile
                                if (hcmref .gt. 0.0 .and. maxval(fafraczint) .le. 0.0) then !revert to prescribed shape profile
                                    call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                        fafraczint)
                                end if
                            else !revert back to using prescribed shape profile
                                call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                    fafraczint)
                            end if
                        end if

! ... calculate leaf area density profile from foliage shape function for output (m2/m3)
                        do k=1, modlays
                            if (zk(k) .gt. 0.0 .and. zk(k) .le. hcmref) then ! above ground level and at/below canopy top
                                if (k .lt. modlays)  then
                                    lad(loc,k) = ((fafraczint(k+1) - fafraczint(k))*lairef)/modres
                                else
                                    lad(loc,k) = lad(loc,modlays-1)
                                end if
                            else
                                lad(loc,k) = 0.0_rk
                            end if
                        end do

! ... calculate zero-plane displacement height/hc and surface (soil+veg) roughness lengths/hc

                        call canopy_zpd(zhc(1:cansublays), fafraczint(1:cansublays), &
                            ubzref, z0ghc, lambdars, cdrag, pai, hcmref, hgtref, &
                            z0ref, vtyperef, lu_opt, z0_opt, d_h(loc), zo_h(loc))

! ... calculate canopy radiation (sunlit and shade) profile

                        call canopy_fsun_clu( fafraczint, lairef, cluref, cszref, fsun)

! ... calculate canopy leaf temperature (sun/shade) profile

                        call canopy_tleaf_lin(zk, hcmref, tmp2mref, fsun, &
                            tleaf_sun, tleaf_shade, tleaf_ave)

! ... calculate canopy Photosynthetic Photon Flux Density (PPFD) (sun/shade) profile

                        call canopy_ppfd_exp(zk, hcmref, vbdsf_averef, vddsf_averef, lairef, fsun, &
                            ppfd_sun, ppfd_shade, ppfd_ave)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                        if (ifcanwind .or. ifcanwaf) then
                            if (rsl_opt .eq. 0) then
                                do k=1, modlays
                                    call canopy_wind_most(hcmref, zk(k), fafraczint(k), ubzref, &
                                        z0ghc, cdrag, pai, hgtref, d_h(loc), zo_h(loc), &
                                        lambdars, canbot(k), cantop(k), canwind(loc,k))
                                end do
                            else
                                write(*,*) 'wrong RSL_OPT namelist option = ', rsl_opt, 'only option = 0 right now'
                                call exit(2)
                            end if

! ... determine midflamepoint and flame height from user or FRP calculation
                            call canopy_flameh(flameh_opt, flameh_set, dx(loc), modres, &
                                frpref, frp_fac, hcmref, lu_opt, vtyperef, flameh_cal, &
                                midflamepoint, flameh(loc))
                            if (firetype .eq. 0) then !forest/sub-canopy firetype
                                if (flameh(loc) .gt. 0.0) then !flameh must be > 0
                                    if (flameh(loc) .le. hcmref) then !only calculate when flameh <= FCH
                                        call canopy_waf(hcmref, lambdars, hgtref, flameh(loc), &
                                            firetype, d_h(loc), zo_h(loc), canbot(midflamepoint), &
                                            cantop(midflamepoint), waf(loc))
                                    else
                                        write(*,*) 'warning...sub-canopy type fire, but flameh > FCH, setting WAF=1'
                                        waf(loc) = 1.0_rk
                                    end if
                                end if
                            else  !grass/crops, above-canopy firetype
                                if (flameh(loc) .gt. 0.0) then !flameh still must be > 0
                                    call canopy_waf(hcmref, lambdars, hgtref, flameh(loc), &
                                        firetype, d_h(loc), zo_h(loc), canbot(midflamepoint), &
                                        cantop(midflamepoint), waf(loc))
                                end if
                            end if
                            if (waf(loc) .gt. 1.0_rk) then !Final check of WAF > 1, must be <=1
                                waf(loc) = 1.0_rk
                            end if
                        end if

! ... user option to calculate in-canopy eddy diffusivities at height z
                        if (ifcaneddy) then
                            do k=1, modlays
                                call canopy_eddyx(hcmref, zk(k), ustref, molref, kz(loc, k))
                            end do
                        end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                        if (ifcanphot) then
                            if (cszref .ge. 0.0_rk) then !only calculate if cell isn't dark
                                call canopy_phot(fafraczint, &
                                    lairef, cluref, cszref, rjcf(loc, :))
                            end if
                        end if

!.......user option to calculate in-canopy leafage influence and assigning LAI as per timestep

                        if (leafage_opt .eq. 0) then
                            !!! Check if the lai_tstep is greater than time_intvl
                            if (lai_tstep .ge. time_intvl) then
                                tsteplai = lai_tstep/86400.0_rk  !convert lai_tstep to time step in # of days needed for bud break in leafage
                                nlaic = nn*time_intvl/lai_tstep  !calculate number of lai tsteps that have elapsed
                                int_nlaic = int(nlaic)
                            else
                                WRITE (*, *) "Error: Input LAI time step cannot be less than model time step...exiting!!!"
                                CALL exit(1)
                            endif
                            ! Initialize pastlai and currentlai based on current timestep
                            if (nn .eq. 1) then
                                currentlai = lairef
                                pastlai = currentlai
                            else if (int_nlaic .gt. int_nlaip) then !only update with each new lai tstep
                                pastlai = currentlai
                                currentlai = lairef
                            endif
                            nlaip = nlaic !set value to compare in next model time step
                            int_nlaip = int(nlaip)
                        end if !leafage_opt = 0 end

!.......user option to calculate historical leaf temperature and PAR for past 24-hours and 240-hours rolling average per timestep
                        !Initialize
                        if (hist_opt .eq. 0) then      !Use instantaneous only
                            ppfd_sun24(loc,:)     = ppfd_sun
                            ppfd_shade24(loc,:)   = ppfd_shade
                            tleaf_sun24(loc,:)    = tleaf_sun
                            tleaf_shade24(loc,:)  = tleaf_shade
                            tleaf_ave24(loc,:)    = tleaf_ave
                            ppfd_sun240(loc,:)    = ppfd_sun
                            ppfd_shade240(loc,:)  = ppfd_shade
                            tleaf_sun240(loc,:)   = tleaf_sun
                            tleaf_shade240(loc,:) = tleaf_shade
                            tleaf_ave240(loc,:)   = tleaf_ave
                            !for AQ, temp, and wind stress factors in biogenics, use instantaneous only (caution!)
                            daily_maxt2m(loc)     = tmp2mref
                            daily_mint2m(loc)     = tmp2mref
                            daily_maxws10m(loc)   = ubzref
                        else if (hist_opt .eq. 1) then  !Try for historical average values
                            ! Calculate weights for running means of historic variables
                            ! DNEWFRAC and DOLDFRAC are the weights given to the current
                            ! and existing value, respectively, when updating running means
                            ! over the last X days. HNEWFRAC and HOLDFRAC are the same but
                            ! for the 24H means.
                            dnewfrac = time_intvl / ( tau_days * 24.0_rk * 3600.0_rk )
                            doldfrac = 1.0_rk - dnewfrac
                            hnewfrac = time_intvl / ( tau_hours * 3600.0_rk )
                            holdfrac = 1.0_rk - hnewfrac

                            if (nn .le. 24) then
                                !Track times for moving time window (running) average below
                                ppfd_sun24_tmp(nn,loc,:)     = ppfd_sun
                                ppfd_shade24_tmp(nn,loc,:)   = ppfd_shade
                                tleaf_sun24_tmp(nn,loc,:)    = tleaf_sun
                                tleaf_shade24_tmp(nn,loc,:)  = tleaf_shade
                                tleaf_ave24_tmp(nn,loc,:)    = tleaf_ave
                                ppfd_sun240_tmp(nn,loc,:)    = ppfd_sun
                                ppfd_shade240_tmp(nn,loc,:)  = ppfd_shade
                                tleaf_sun240_tmp(nn,loc,:)   = tleaf_sun
                                tleaf_shade240_tmp(nn,loc,:) = tleaf_shade
                                tleaf_ave240_tmp(nn,loc,:)   = tleaf_ave
                                tmp2mref_tmp(nn,loc)         = tmp2mref
                                ubzref_tmp(nn,loc)           = ubzref

                                !TODO:  Restart capability needed to get past leaf temp and PAR if avaialble
                                !For now, if <= 24 hours then use current instantaneous
                                ppfd_sun24(loc,:)     = ppfd_sun
                                ppfd_shade24(loc,:)   = ppfd_shade
                                tleaf_sun24(loc,:)    = tleaf_sun
                                tleaf_shade24(loc,:)  = tleaf_shade
                                tleaf_ave24(loc,:)    = tleaf_ave
                                ppfd_sun240(loc,:)    = ppfd_sun
                                ppfd_shade240(loc,:)  = ppfd_shade
                                tleaf_sun240(loc,:)   = tleaf_sun
                                tleaf_shade240(loc,:) = tleaf_shade
                                tleaf_ave240(loc,:)   = tleaf_ave
                                !for AQ, temp, and wind stress factors in biogenics, still use instantaneous only (caution!)
                                daily_maxt2m(loc)     = tmp2mref
                                daily_mint2m(loc)     = tmp2mref
                                daily_maxws10m(loc)   = ubzref

                            else if (nn .le. 240) then
                                !Drop nth value into end time slice
                                ppfd_sun24_tmp(25,loc,:)     = ppfd_sun
                                ppfd_shade24_tmp(25,loc,:)   = ppfd_shade
                                tleaf_sun24_tmp(25,loc,:)    = tleaf_sun
                                tleaf_shade24_tmp(25,loc,:)  = tleaf_shade
                                tleaf_ave24_tmp(25,loc,:)    = tleaf_ave
                                ppfd_sun240_tmp(nn,loc,:)    = ppfd_sun
                                ppfd_shade240_tmp(nn,loc,:)  = ppfd_shade
                                tleaf_sun240_tmp(nn,loc,:)   = tleaf_sun
                                tleaf_shade240_tmp(nn,loc,:) = tleaf_shade
                                tleaf_ave240_tmp(nn,loc,:)   = tleaf_ave
                                tmp2mref_tmp(25,loc)         = tmp2mref
                                ubzref_tmp(25,loc)           = ubzref

                                !Calculate 24 and 240 hour averages
                                ppfd_sun24(loc,:)     = sum(ppfd_sun24_tmp(1:24,loc,:),1)/24.0
                                ppfd_shade24(loc,:)   = sum(ppfd_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_sun24(loc,:)    = sum(tleaf_sun24_tmp(1:24,loc,:),1)/24.0
                                tleaf_shade24(loc,:)  = sum(tleaf_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_ave24(loc,:)    = sum(tleaf_ave24_tmp(1:24,loc,:),1)/24.0
                                !Use 24 hour restarts until 240 restarts become available
                                ppfd_sun240(loc,:)    = sum(ppfd_sun24_tmp(1:24,loc,:),1)/24.0
                                ppfd_shade240(loc,:)  = sum(ppfd_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_sun240(loc,:)   = sum(tleaf_sun24_tmp(1:24,loc,:),1)/24.0
                                tleaf_shade240(loc,:) = sum(tleaf_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_ave240(loc,:)   = sum(tleaf_ave24_tmp(1:24,loc,:),1)/24.0
                                !Updated rolling 24 hour (hourly, short term) and 240 hour (daily, long-term) averages
                                ppfd_sun24(loc,:)     = ( holdfrac * ppfd_sun24(loc,:) )     + ( hnewfrac * ppfd_sun )
                                ppfd_shade24(loc,:)   = ( holdfrac * ppfd_shade24(loc,:) )   + ( hnewfrac * ppfd_shade )
                                tleaf_sun24(loc,:)    = ( holdfrac * tleaf_sun24(loc,:) )    + ( hnewfrac * tleaf_sun )
                                tleaf_shade24(loc,:)  = ( holdfrac * tleaf_shade24(loc,:) )  + ( hnewfrac * tleaf_shade )
                                tleaf_ave24(loc,:)    = ( holdfrac * tleaf_ave24(loc,:) )    + ( hnewfrac * tleaf_ave )
                                ppfd_sun240(loc,:)     = ( holdfrac * ppfd_sun24(loc,:) )     + ( hnewfrac * ppfd_sun )
                                ppfd_shade240(loc,:)   = ( holdfrac * ppfd_shade24(loc,:) )   + ( hnewfrac * ppfd_shade )
                                tleaf_sun240(loc,:)    = ( holdfrac * tleaf_sun24(loc,:) )    + ( hnewfrac * tleaf_sun )
                                tleaf_shade240(loc,:)  = ( holdfrac * tleaf_shade24(loc,:) )  + ( hnewfrac * tleaf_shade )
                                tleaf_ave240(loc,:)    = ( holdfrac * tleaf_ave24(loc,:) )    + ( hnewfrac * tleaf_ave )
                                !Take daily max/min for temperature and daily max for wind speed
                                daily_maxt2m(loc)        = maxval(tmp2mref_tmp(1:24,loc))
                                daily_mint2m(loc)        = minval(tmp2mref_tmp(1:24,loc))
                                daily_maxws10m(loc)      = maxval(ubzref_tmp(1:24,loc))

                                !Shift time window to make room for new nth value
                                ppfd_sun24_tmp(1:24,loc,:)     = ppfd_sun24_tmp(2:25,loc,:)
                                ppfd_shade24_tmp(1:24,loc,:)   = ppfd_shade24_tmp(2:25,loc,:)
                                tleaf_sun24_tmp(1:24,loc,:)    = tleaf_sun24_tmp(2:25,loc,:)
                                tleaf_shade24_tmp(1:24,loc,:)  = tleaf_shade24_tmp(2:25,loc,:)
                                tleaf_ave24_tmp(1:24,loc,:)    = tleaf_ave24_tmp(2:25,loc,:)
                                tmp2mref_tmp(1:24,loc)         = tmp2mref_tmp(2:25,loc)
                                ubzref_tmp(1:24,loc)           = ubzref_tmp(2:25,loc)

                            else
                                !Drop nth value into end time slice
                                ppfd_sun24_tmp(25,loc,:)     = ppfd_sun
                                ppfd_shade24_tmp(25,loc,:)   = ppfd_shade
                                tleaf_sun24_tmp(25,loc,:)    = tleaf_sun
                                tleaf_shade24_tmp(25,loc,:)  = tleaf_shade
                                tleaf_ave24_tmp(25,loc,:)    = tleaf_ave
                                ppfd_sun240_tmp(241,loc,:)    = ppfd_sun
                                ppfd_shade240_tmp(241,loc,:)  = ppfd_shade
                                tleaf_sun240_tmp(241,loc,:)   = tleaf_sun
                                tleaf_shade240_tmp(241,loc,:) = tleaf_shade
                                tleaf_ave240_tmp(241,loc,:)   = tleaf_ave
                                tmp2mref_tmp(25,loc)         = tmp2mref
                                ubzref_tmp(25,loc)           = ubzref

                                !Calculate 24 and 240 hour averages
                                ppfd_sun24(loc,:)     = sum(ppfd_sun24_tmp(1:24,loc,:),1)/24.0
                                ppfd_shade24(loc,:)   = sum(ppfd_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_sun24(loc,:)    = sum(tleaf_sun24_tmp(1:24,loc,:),1)/24.0
                                tleaf_shade24(loc,:)  = sum(tleaf_shade24_tmp(1:24,loc,:),1)/24.0
                                tleaf_ave24(loc,:)    = sum(tleaf_ave24_tmp(1:24,loc,:),1)/24.0
                                ppfd_sun240(loc,:)    = sum(ppfd_sun240_tmp(1:240,loc,:),1)/240.0
                                ppfd_shade240(loc,:)  = sum(ppfd_shade240_tmp(1:240,loc,:),1)/240.0
                                tleaf_sun240(loc,:)   = sum(tleaf_sun240_tmp(1:240,loc,:),1)/240.0
                                tleaf_shade240(loc,:) = sum(tleaf_shade240_tmp(1:240,loc,:),1)/240.0
                                tleaf_ave240(loc,:)   = sum(tleaf_ave240_tmp(1:240,loc,:),1)/240.0
                                !Updated rolling 24 hour (hourly, short term) and 240 hour (daily, long-term) averages
                                ppfd_sun24(loc,:)     = ( holdfrac * ppfd_sun24(loc,:) )     + ( hnewfrac * ppfd_sun )
                                ppfd_shade24(loc,:)   = ( holdfrac * ppfd_shade24(loc,:) )   + ( hnewfrac * ppfd_shade )
                                tleaf_sun24(loc,:)    = ( holdfrac * tleaf_sun24(loc,:) )    + ( hnewfrac * tleaf_sun )
                                tleaf_shade24(loc,:)  = ( holdfrac * tleaf_shade24(loc,:) )  + ( hnewfrac * tleaf_shade )
                                tleaf_ave24(loc,:)    = ( holdfrac * tleaf_ave24(loc,:) )    + ( hnewfrac * tleaf_ave )
                                ppfd_sun240(loc,:)    = ( doldfrac * ppfd_sun240(loc,:) )    + ( dnewfrac * ppfd_sun )
                                ppfd_shade240(loc,:)  = ( doldfrac * ppfd_shade240(loc,:) )  + ( dnewfrac * ppfd_shade )
                                tleaf_sun240(loc,:)   = ( doldfrac * tleaf_sun240(loc,:) )   + ( dnewfrac * tleaf_sun )
                                tleaf_shade240(loc,:) = ( doldfrac * tleaf_shade240(loc,:) ) + ( dnewfrac * tleaf_shade )
                                tleaf_ave240(loc,:)   = ( doldfrac * tleaf_ave240(loc,:) )   + ( dnewfrac * tleaf_ave )
                                !Take daily max/min for temperature and daily max for wind speed
                                daily_maxt2m(loc)        = maxval(tmp2mref_tmp(1:24,loc))
                                daily_mint2m(loc)        = minval(tmp2mref_tmp(1:24,loc))
                                daily_maxws10m(loc)      = maxval(ubzref_tmp(1:24,loc))

                                !Shift time window to make room for new nth value
                                ppfd_sun24_tmp(1:24,loc,:)     = ppfd_sun24_tmp(2:25,loc,:)
                                ppfd_shade24_tmp(1:24,loc,:)   = ppfd_shade24_tmp(2:25,loc,:)
                                tleaf_sun24_tmp(1:24,loc,:)    = tleaf_sun24_tmp(2:25,loc,:)
                                tleaf_shade24_tmp(1:24,loc,:)  = tleaf_shade24_tmp(2:25,loc,:)
                                tleaf_ave24_tmp(1:24,loc,:)    = tleaf_ave24_tmp(2:25,loc,:)
                                ppfd_sun240_tmp(1:240,loc,:)    = ppfd_sun240_tmp(2:241,loc,:)
                                ppfd_shade240_tmp(1:240,loc,:)  = ppfd_shade240_tmp(2:241,loc,:)
                                tleaf_sun240_tmp(1:240,loc,:)   = tleaf_sun240_tmp(2:241,loc,:)
                                tleaf_shade240_tmp(1:240,loc,:) = tleaf_shade240_tmp(2:241,loc,:)
                                tleaf_ave240_tmp(1:240,loc,:)   = tleaf_ave240_tmp(2:241,loc,:)
                                tmp2mref_tmp(1:24,loc)          = tmp2mref_tmp(2:25,loc)
                                ubzref_tmp(1:24,loc)            = ubzref_tmp(2:25,loc)
                            end if
                        else
                            write(*,*) 'wrong HIST_OPT namelist option = ', hist_opt, 'only option = 0 or 1'
                            call exit(2)
                        end if
! ... user option to calculate in-canopy biogenic emissions
                        if (ifcanbio) then
                            if (cluref .gt. 0.0_rk) then
                                !ISOP
                                if (biospec_opt == 0 .or. biospec_opt == 1) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 1, emi_isop(loc,:))
                                end if
                                !MYRC
                                if (biospec_opt == 0 .or. biospec_opt == 2) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 2, emi_myrc(loc,:))
                                end if
                                !SABI
                                if (biospec_opt == 0 .or. biospec_opt == 3) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 3, emi_sabi(loc,:))
                                end if
                                !LIMO
                                if (biospec_opt == 0 .or. biospec_opt == 4) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 4, emi_limo(loc,:))
                                end if
                                !CARE
                                if (biospec_opt == 0 .or. biospec_opt == 5) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 5, emi_care(loc,:))
                                end if
                                !OCIM
                                if (biospec_opt == 0 .or. biospec_opt == 6) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 6, emi_ocim(loc,:))
                                end if
                                !BPIN
                                if (biospec_opt == 0 .or. biospec_opt == 7) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 7, emi_bpin(loc,:))
                                end if
                                !APIN
                                if (biospec_opt == 0 .or. biospec_opt == 8) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 8, emi_apin(loc,:))
                                end if
                                !MONO
                                if (biospec_opt == 0 .or. biospec_opt == 9) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 9, emi_mono(loc,:))
                                end if
                                !FARN
                                if (biospec_opt == 0 .or. biospec_opt == 10) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 10, emi_farn(loc,:))
                                end if
                                !CARY
                                if (biospec_opt == 0 .or. biospec_opt == 11) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 11, emi_cary(loc,:))
                                end if
                                !SESQ
                                if (biospec_opt == 0 .or. biospec_opt == 12) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 12, emi_sesq(loc,:))
                                end if
                                !MBOL
                                if (biospec_opt == 0 .or. biospec_opt == 13) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 13, emi_mbol(loc,:))
                                end if
                                !METH
                                if (biospec_opt == 0 .or. biospec_opt == 14) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 14, emi_meth(loc,:))
                                end if
                                !ACET
                                if (biospec_opt == 0 .or. biospec_opt == 15) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 15, emi_acet(loc,:))
                                end if
                                !CO
                                if (biospec_opt == 0 .or. biospec_opt == 16) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 16, emi_co(loc,:))
                                end if
                                !BIDI VOC
                                if (biospec_opt == 0 .or. biospec_opt == 17) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 17, emi_bvoc(loc,:))
                                end if
                                !Stress VOC
                                if (biospec_opt == 0 .or. biospec_opt == 18) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 18, emi_svoc(loc,:))
                                end if
                                !Other VOC
                                if (biospec_opt == 0 .or. biospec_opt == 19) then
                                    call canopy_bio(zk, fafraczint, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade,&
                                        ppfd_sun24(loc,:), ppfd_shade24(loc,:), &
                                        tleaf_ave24(loc,:), ppfd_sun240(loc,:), ppfd_shade240(loc,:), &
                                        tleaf_ave240(loc,:), tka(loc,:), dswrfref, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        loss_opt, loss_set, loss_ind, lifetime, ustref, &
                                        soim_opt, soilw1ref, soilw2ref, soilw3ref, soilw4ref, &
                                        soild1, soild2, soild3, soild4, wiltref, aq_opt, w126_set, ozone_w126ref, &
                                        ht_opt, lt_opt, hw_opt, daily_maxt2m(loc), daily_mint2m(loc), &
                                        daily_maxws10m(loc), &
                                        modlays, 19, emi_ovoc(loc,:))
                                end if
                            else
                                if (biospec_opt == 0 .or. biospec_opt == 1) then
                                    emi_isop(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 2) then
                                    emi_myrc(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 3) then
                                    emi_sabi(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 4) then
                                    emi_limo(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 5) then
                                    emi_care(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 6) then
                                    emi_ocim(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 7) then
                                    emi_bpin(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 8) then
                                    emi_apin(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 9) then
                                    emi_mono(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 10) then
                                    emi_farn(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 11) then
                                    emi_cary(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 12) then
                                    emi_sesq(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 13) then
                                    emi_mbol(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 14) then
                                    emi_meth(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 15) then
                                    emi_acet(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 16) then
                                    emi_co(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 17) then
                                    emi_bvoc(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 18) then
                                    emi_svoc(loc,:) = 0.0_rk
                                end if
                                if (biospec_opt == 0 .or. biospec_opt == 19) then
                                    emi_ovoc(loc,:) = 0.0_rk
                                end if
                            end if
                        end if

! ... user option to calculate in-canopy dry deposition velocity
                        if (ifcanddepgas ) then
                            if (ifcanwind) then   !ubar needed for rbl
                                if (chemmechgas_opt .eq. 0) then   !RACM2
                                    if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                        !At this point must have canopy height set by thresholds, so Ra calculation is valid
                                        !Calculate Bulk Richardson Number for stability dependent aerodynamic resistance
                                        rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                            d_h(loc)*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk)
                                        ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, d_h(loc)*hcmref*100.0_rk, &
                                            hcmref*100.0_rk, rib)
                                        ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum,and assume constant Ra
                                        !downward at all levels through canopy to ground
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 1, ddep_no(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra, 1, ddep_no(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 1,ddep_no(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 2, ddep_no2(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  2, ddep_no2(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 2,ddep_no2(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 3, ddep_o3(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  3, ddep_o3(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 3,ddep_o3(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 4, ddep_hono(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  4, ddep_hono(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 4,ddep_hono(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 5, ddep_hno4(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  5, ddep_hno4(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 5,ddep_hno4(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 6, ddep_hno3(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  6, ddep_hno3(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 6,ddep_hno3(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 7, ddep_n2o5(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  7, ddep_n2o5(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 7,ddep_n2o5(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 8, ddep_co(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  8, ddep_co(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 8,ddep_co(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 9, ddep_h2o2(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  9, ddep_h2o2(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 9,ddep_h2o2(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref,  ra,10, ddep_ch4(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  10, ddep_ch4(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 10,ddep_ch4(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 11, ddep_mo2(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  11, ddep_mo2(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 11,ddep_mo2(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 12, ddep_op1(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref,  ra, 12, ddep_op1(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 12,ddep_op1(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 13, ddep_moh(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  13, ddep_moh(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 13,ddep_moh(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 14, ddep_no3(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  14, ddep_no3(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 14,ddep_no3(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt ==15) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 15, ddep_o3p(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  15, ddep_o3p(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 15,ddep_o3p(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 16, ddep_o1d(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  16, ddep_o1d(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 16,ddep_o1d(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 17, ddep_ho(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  17, ddep_ho(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 17,ddep_ho(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 18, ddep_ho2(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  18, ddep_ho2(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 18,ddep_ho2(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 19, ddep_ora1(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  19, ddep_ora1(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 19,ddep_ora1(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 20, ddep_hac(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  20, ddep_hac(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 20,ddep_hac(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 21, ddep_paa(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  21, ddep_paa(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 21,ddep_paa(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 22, ddep_dhmob(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  22, ddep_dhmob(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 22,ddep_dhmob(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 23, ddep_hpald(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  23, ddep_hpald(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 23,ddep_hpald(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 24, ddep_ishp(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  24, ddep_ishp(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 24,ddep_ishp(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 25, ddep_iepox(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  25, ddep_iepox(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 25,ddep_iepox(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 26, ddep_propnn(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  26, ddep_propnn(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 26,ddep_propnn(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 27, ddep_isopnb(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  27, ddep_isopnb(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 27,ddep_isopnb(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 28, ddep_isopnd(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  28, ddep_isopnd(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 28,ddep_isopnd(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 29, ddep_macrn(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  29, ddep_macrn(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 29,ddep_macrn(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 30, ddep_mvkn(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  30, ddep_mvkn(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 30,ddep_mvkn(loc,1))
                                            end if
                                        endif
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                            call canopy_gas_drydep_zhang(chemmechgas_opt,chemmechgas_tot, &
                                                zk, hcmref, tka(loc,:), pressa(loc,:), &
                                                relhuma(loc,:), fsun, ppfd_sun, ppfd_shade, canwind(loc,:),  &
                                                dswrfref, ra, 31, ddep_isnp(loc,:))   ! [cm/s]
                                        endif
                                        !soil gas dry depostion at level 1, i.e., z=0
                                        if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                            if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                                call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                                    soilt1ref, pressa(loc,1), canwind(loc,2),  &
                                                    soilcat_opt,sotypref,soild1, soilw1ref, ra,  31, ddep_isnp(loc,1))   ! [cm/s]
                                            else !snow covered  or ice covered
                                                call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                                    canwind(loc,2),ra, 31,ddep_isnp(loc,1))
                                            end if
                                        endif
                                    else
                                        write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                        write(*,*) ' in namelist...exiting'
                                        write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                        call exit(2)
                                    end if
                                else
                                    write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                                    write(*,*)  'Set chemmechgas_opt to only 0 (RACM2) for now'
                                    call exit(2)
                                end if
                            else
                                write(*,*)  'Wrong IfCanWind choice of ', ifcanwind, ' in namelist...exiting'
                                write(*,*)  'Set IfCanwind to True to use IfCanDDepGas'
                                call exit(2)
                            end if
                        end if

                        ! --- Sub-canopy aerosol dry deposition through vegetative canopies (Katul et al. 2010) ---
                        if (ifcanaeroddep) then
                            if (ifcanwind) then !ubar needed for rbl
                                rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                                    d_h(loc)*hcmref, (hcmref+hgtref))
                                ra=rav(ubzref, (hcmref+hgtref), d_h(loc)*hcmref, &
                                    hcmref, rib)
                                ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                                call canopy_aero_ddep_subveg(modlays, zk, hcmref, lad(loc,:), canwind(loc,:), &
                                    aeroddep_diam, aeroddep_rho, tka(loc,:), pressa(loc,:), aeroddep_opt, &
                                    ra, modres, ustref, vtyperef, vdep_aero(loc,:)) ! [cm/s]
                            else
                                write(*,*)  'Wrong IfCanWind choice of ', ifcanwind, ' in namelist...exiting'
                                write(*,*)  'Set IfCanwind to True to use IfCanAeroDDep'
                                call exit(2)
                            end if
                            !  --- Sub-canopy aerosol dry deposition to soil (i.e., first model layer) beneath canopy (Pleim et al. 2022) ---
                            rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                                d_h(loc)*hcmref, (hcmref+hgtref))
                            ra=rav(ubzref, (hcmref+hgtref), d_h(loc)*hcmref, &
                                hcmref, rib)
                            ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                            call canopy_aero_ddep_pleim2022(ustref, ra, canwind(loc,2), &
                                aeroddep_diam, aeroddep_rho, tka(loc,1), pressa(loc,1), 2, vdep_aero(loc,1))
                        end if

                    else
                        if (biospec_opt == 0 .or. biospec_opt == 1) then
                            emi_isop(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 2) then
                            emi_myrc(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 3) then
                            emi_sabi(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 4) then
                            emi_limo(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 5) then
                            emi_care(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 6) then
                            emi_ocim(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 7) then
                            emi_bpin(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 8) then
                            emi_apin(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 9) then
                            emi_mono(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 10) then
                            emi_farn(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 11) then
                            emi_cary(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 12) then
                            emi_sesq(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 13) then
                            emi_mbol(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 14) then
                            emi_meth(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 15) then
                            emi_acet(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 16) then
                            emi_co(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 17) then
                            emi_bvoc(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 18) then
                            emi_svoc(loc,:) = 0.0_rk
                        end if
                        if (biospec_opt == 0 .or. biospec_opt == 19) then
                            emi_ovoc(loc,:) = 0.0_rk
                        end if
                    end if !Contiguous Canopy

                else if (vtyperef .eq. 15 .or. vtyperef .eq. 16 .or. vtyperef .eq. 20) then !Barren/Sparsely Vegetated  or Barren
                    !Tundra or Snow/Ice
! ... user option to calculate dry deposition velocity...for land use outside of vegetated canopies
                    if (ifcanddepgas ) then
                        if (chemmechgas_opt .eq. 0)  then   !RACM2
                            if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (hcmref .le. 0.0) then
                                    hcmref = 0.01!m   (set to just above ground @ 1 cm; leads to z0 ~ 0.001 m
                                endif                ! needed for Ra calculation)
                                rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                    0.75_rk*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk) !assume d~3/4*Hc
                                ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, 0.75_rk*hcmref*100.0_rk, &
                                    hcmref*100.0_rk, rib)
                                ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  1, ddep_no(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 1,ddep_no(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref,  ra, 2, ddep_no2(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 2,ddep_no2(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  3, ddep_o3(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 3,ddep_o3(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  4, ddep_hono(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 4,ddep_hono(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  5, ddep_hno4(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 5,ddep_hno4(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  6, ddep_hno3(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 6,ddep_hno3(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  7, ddep_n2o5(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 7,ddep_n2o5(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  8, ddep_co(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 8,ddep_co(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  9, ddep_h2o2(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 9,ddep_h2o2(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  10, ddep_ch4(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 10,ddep_ch4(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  11, ddep_mo2(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 11,ddep_mo2(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref,  ra, 12, ddep_op1(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 12,ddep_op1(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  13, ddep_moh(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 13,ddep_moh(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  14, ddep_no3(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 14,ddep_no3(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  15, ddep_o3p(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 15,ddep_o3p(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  16, ddep_o1d(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 16,ddep_o1d(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  17, ddep_ho(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 17,ddep_ho(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  18, ddep_ho2(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 18,ddep_ho2(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  19, ddep_ora1(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 19,ddep_ora1(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  20, ddep_hac(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 20,ddep_hac(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  21, ddep_paa(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 21,ddep_paa(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  22, ddep_dhmob(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 22,ddep_dhmob(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  23, ddep_hpald(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 23,ddep_hpald(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  24, ddep_ishp(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 24,ddep_ishp(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  25, ddep_iepox(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 25,ddep_iepox(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  26, ddep_propnn(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 26,ddep_propnn(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  27, ddep_isopnb(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 27,ddep_isopnb(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  28, ddep_isopnd(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 28,ddep_isopnd(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  29, ddep_macrn(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 29,ddep_macrn(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  30, ddep_mvkn(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 30,ddep_mvkn(loc,1))
                                    end if
                                endif
                                !soil (barren land) gas dry depostion at level 1, i.e., z=0
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_soil(chemmechgas_opt,chemmechgas_tot, &
                                            soilt1ref, pressa(loc,1), ubzref,  &
                                            soilcat_opt,sotypref,soild1, soilw1ref, ra,  31, ddep_isnp(loc,1))   ! [cm/s]
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 31,ddep_isnp(loc,1))
                                    end if
                                endif
                            else
                                write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                write(*,*)  ' in namelist...exiting'
                                write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                call exit(2)
                            end if
                        else
                            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                            write(*,*)  'Set chemmechgas_opt = 0 (RACM2) for now'
                            call exit(2)
                        end if
                    end if
                    ! --- Sub-canopy aerosol dry deposition over barren/soil surfaces (Pleim et al. 2022) ---
                    if (ifcanaeroddep) then
                        rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                            d_h(loc)*hcmref, (hcmref+hgtref))
                        ra=rav(ubzref, (hcmref+hgtref), d_h(loc)*hcmref, &
                            hcmref, rib)
                        ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                        call canopy_aero_ddep_pleim2022(ustref, ra, ubzref, &
                            aeroddep_diam, aeroddep_rho, tmp2mref, pressfcref, 2, vdep_aero(loc,1))
                    end if
                else if (vtyperef .eq. 13) then !Urban and Built Up
! ... user option to calculate dry deposition velocity...for land use outside of vegetated canopies
                    if (ifcanddepgas ) then
                        if (chemmechgas_opt .eq. 0)  then   !RACM2
                            if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                !urban gas dry depostion at level 1, i.e., z=0
                                if (hcmref .le. 0.0) then
                                    hcmref = 0.01!m   (set to just above ground @ 1 cm; leads to z0 ~ 0.001 m
                                endif                ! needed for Ra calculation)
                                rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                    0.75_rk*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk) !assume d~3/4*Hc
                                ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, 0.75_rk*hcmref*100.0_rk, &
                                    hcmref*100.0_rk, rib)
                                ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra, 1,ddep_no(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 1,ddep_no(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra, 2,ddep_no2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 2,ddep_no2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra, 3,ddep_o3(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 3,ddep_o3(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,4,ddep_hono(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 4,ddep_hono(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,5,ddep_hno4(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 5,ddep_hno4(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,6,ddep_hno3(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 6,ddep_hno3(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,7,ddep_n2o5(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 7,ddep_n2o5(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,8,ddep_co(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 8,ddep_co(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,9,ddep_h2o2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 9,ddep_h2o2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,10,ddep_ch4(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 10,ddep_ch4(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,11,ddep_mo2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 11,ddep_mo2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,12,ddep_op1(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 12,ddep_op1(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,13,ddep_moh(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 13,ddep_moh(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,14,ddep_no3(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 14,ddep_no3(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,15,ddep_o3p(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 15,ddep_o3p(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,16,ddep_o1d(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 16,ddep_o1d(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,17,ddep_ho(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 17,ddep_ho(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,18,ddep_ho2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 18,ddep_ho2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,19,ddep_ora1(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 19,ddep_ora1(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,20,ddep_hac(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 20,ddep_hac(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,21,ddep_paa(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 21,ddep_paa(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,22,ddep_dhmob(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 22,ddep_dhmob(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,23,ddep_hpald(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 23,ddep_hpald(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,24,ddep_ishp(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 24,ddep_ishp(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,25,ddep_iepox(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 25,ddep_iepox(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,26,ddep_propnn(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 26,ddep_propnn(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,27,ddep_isopnb(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 27,ddep_isopnb(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,28,ddep_isopnd(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 28,ddep_isopnd(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,29,ddep_macrn(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 29,ddep_macrn(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,30,ddep_mvkn(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 30,ddep_mvkn(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_urban(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref,tmp2mref,gamma_set, ra,31,ddep_isnp(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 31,ddep_isnp(loc,1))
                                    end if
                                endif
                            else
                                write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                write(*,*)  ' in namelist...exiting'
                                write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                call exit(2)
                            end if
                        else
                            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                            write(*,*)  'Set chemmechgas_opt = 0 (RACM2) for now'
                            call exit(2)
                        end if
                    end if
                    ! --- Sub-canopy aerosol dry deposition over urban surfaces (Pleim et al. 2022) ---
                    if (ifcanaeroddep) then
                        rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                            d_h(loc)*hcmref, (hcmref+hgtref))
                        ra=rav(ubzref, (hcmref+hgtref), d_h(loc)*hcmref, &
                            hcmref, rib)
                        ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                        call canopy_aero_ddep_pleim2022(ustref, ra, ubzref, &
                            aeroddep_diam, aeroddep_rho, tmp2mref, pressfcref, 1, vdep_aero(loc,1))
                    end if
                else if (vtyperef .eq. 0) then !Water from FV3 (usually vtype = 17 for water)
! ... user option to calculate dry deposition velocity...for land use outside of vegetated canopies
                    if (ifcanddepgas ) then
                        if (chemmechgas_opt .eq. 0)  then   !RACM2
                            if (chemmechgas_tot .eq. 31) then   !RACM2=31 total gas species including transport
                                !water gas dry depostion at level 1, i.e., z=0
                                if (hcmref .le. 0.0) then
                                    hcmref = 0.01!m   (set to just above ground @ 1 cm; leads to z0 ~ 0.001 m
                                endif                ! needed for Ra calculation)
                                rib=calcrib(tmp2mref, tmpsfcref, ubzref*100.0_rk, &
                                    0.75_rk*hcmref*100.0_rk, (hcmref+hgtref)*100.0_rk) !assume d~3/4*Hc
                                ra=rav(ubzref*100.0_rk, (hcmref+hgtref)*100.0_rk, 0.75_rk*hcmref*100.0_rk, &
                                    hcmref*100.0_rk, rib)
                                ra = max(ramin_set/100.0_rk,ra) !Bound Ra minimum
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,1,ddep_no(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 1,ddep_no(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,2,ddep_no2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 2,ddep_no2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,3,ddep_o3(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 3,ddep_o3(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,4,ddep_hono(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 4,ddep_hono(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,5,ddep_hno4(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 5,ddep_hno4(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,6,ddep_hno3(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 6,ddep_hno3(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,7,ddep_n2o5(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 7,ddep_n2o5(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,8,ddep_co(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 8,ddep_co(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,9,ddep_h2o2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 9,ddep_h2o2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,10,ddep_ch4(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 10,ddep_ch4(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,11,ddep_mo2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 11,ddep_mo2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,12,ddep_op1(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 12,ddep_op1(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,13,ddep_moh(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 13,ddep_moh(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,14,ddep_no3(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 14,ddep_no3(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,15,ddep_o3p(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 15,ddep_o3p(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,16,ddep_o1d(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 16,ddep_o1d(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,17,ddep_ho(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 17,ddep_ho(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,18,ddep_ho2(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 18,ddep_ho2(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,19,ddep_ora1(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 19,ddep_ora1(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,20,ddep_hac(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 20,ddep_hac(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,21,ddep_paa(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 21,ddep_paa(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,22,ddep_dhmob(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 22,ddep_dhmob(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,23,ddep_hpald(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 23,ddep_hpald(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,24,ddep_ishp(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 24,ddep_ishp(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,25,ddep_iepox(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 25,ddep_iepox(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,26,ddep_propnn(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 26,ddep_propnn(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,27,ddep_isopnb(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 27,ddep_isopnb(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,28,ddep_isopnd(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 28,ddep_isopnd(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,29,ddep_macrn(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 29,ddep_macrn(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,30,ddep_mvkn(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 30,ddep_mvkn(loc,1))
                                    end if
                                endif
                                if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                                    if (snowc_averef .le. snowc_set .and. icec_averef*100.0_rk .le. icec_set) then
                                        call canopy_gas_drydep_water(chemmechgas_opt,chemmechgas_tot, &
                                            tmp2mref,spfh2mref,ustref, ra,31,ddep_isnp(loc,1))
                                    else !snow covered  or ice covered
                                        call canopy_gas_drydep_snow(chemmechgas_opt,chemmechgas_tot, &
                                            ubzref, ra, 31,ddep_isnp(loc,1))
                                    end if
                                endif
                            else
                                write(*,*)  'Wrong number of chemical species of ', chemmechgas_tot
                                write(*,*)  ' in namelist...exiting'
                                write(*,*)  'Set chemmechgas_tot = 31 for RACM2'
                                call exit(2)
                            end if
                        else
                            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
                            write(*,*)  'Set chemmechgas_opt = 0 (RACM2) for now'
                            call exit(2)
                        end if
                    end if
                    ! --- Sub-canopy aerosol dry deposition over default surfaces including water (Pleim et al. 2022) ---
                    if (ifcanaeroddep) then
                        rib=calcrib(tmp2mref, tmpsfcref, ubzref, &
                            d_h(loc)*hcmref, (hcmref+hgtref))
                        ra=rav(ubzref, (hcmref+hgtref), d_h(loc)*hcmref, &
                            hcmref, rib)
                        ra = max(ramin_set,ra) !Bound Ra minimum,and assume constant Ra
                        call canopy_aero_ddep_pleim2022(ustref, ra, ubzref, &
                            aeroddep_diam, aeroddep_rho, tmp2mref, pressfcref, 3, vdep_aero(loc,1))
                    end if
                else
                    write(*,*)  'Warning VIIRS/MODIS VTYPE ', vtyperef, ' is not supported...continue'
                end if   !Vegetation types
            else
                write(*,*)  'Wrong LU_OPT choice of ', lu_opt, ' in namelist...exiting'
                call exit(2)

            end if       !Landuse Options
        end do  !Lat*Lon

    else
        write(*,*)  'Wrong INFMT_OPT choice of ', infmt_opt, ' in namelist...exiting'
        call exit(2)

    end if !Input Format (1D or 2D text or NetCDF)

END SUBROUTINE canopy_calcs
```


