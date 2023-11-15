
SUBROUTINE canopy_calcs(nn)

!-------------------------------------------------------------------------------
! Name:     Contains main canopy calculations dependent on canopy conditions
! Purpose:  Contains main canopy calculations dependent on canopy conditions
! Revised:  06 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    use canopy_const_mod, ONLY: rk      !constants for canopy models
    use canopy_coord_mod      !main canopy coordinate descriptions
    use canopy_canopts_mod    !main canopy option descriptions
    use canopy_canmet_mod     !main canopy met/sfc input descriptions
    use canopy_canvars_mod    !main canopy variables descriptions
    use canopy_utils_mod      !main canopy utilities
    use canopy_dxcalc_mod     !main canopy dx calculation
    use canopy_profile_mod    !main canopy foliage profile routines
    use canopy_var3din_mod    !main canopy 3D variable in routines
    use canopy_rad_mod        !main canopy radiation sunlit/shaded routines
    use canopy_tleaf_mod      !main canopy leaf temperature sunlit/shaded routines
    use canopy_wind_mod       !main canopy components
    use canopy_fire_mod
    use canopy_phot_mod
    use canopy_eddy_mod
    use canopy_bioemi_mod

    IMPLICIT NONE

    INTEGER,     INTENT( IN )       :: nn         ! Input time step

    !Local variables
    integer i,j,k,loc
    ! LAI variables for Leaf Age factor calculations
    REAL(rk) :: pastlai       !Past LAI [cm2/cm2]
    REAL(rk), save :: currentlai    ! Current LAI [cm2/cm2]  (saved from one timestep to the next)
    REAL(rk) :: tsteplai  !Number of days between the past and current LAI
    !REAL(rk) :: tabovecanopy  ! Above Canopy Temp (assigned = tmp2mref ), done in canopy_bioemi_mod.F90




    write(*,*)  'Calculating Canopy Parameters'
    write(*,*)  '-------------------------------'

    if (infmt_opt .eq. 0) then !Main input format is 2D NetCDF and output will be 2D NetCDf

        if (ifcanwind .or. ifcanwaf) then !only calculate if canopy wind or WAF option
            call canopy_calcdx_2d(dx_opt, dx_set, nlat, nlon, variables_2d%lat, &
                variables_2d%lon, dx_2d)
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

                hcmref       = variables_2d(i,j)%fh
                uref         = variables_2d(i,j)%ugrd10m
                vref         = variables_2d(i,j)%vgrd10m
                cluref       = variables_2d(i,j)%clu
                lairef       = variables_2d(i,j)%lai
                vtyperef     = variables_2d(i,j)%vtype
                ffracref     = variables_2d(i,j)%ffrac
                ustref       = variables_2d(i,j)%fricv
                cszref       = variables_2d(i,j)%csz
                z0ref        = variables_2d(i,j)%sfcr
                molref       = variables_2d(i,j)%mol
                frpref       = variables_2d(i,j)%frp
                hgtref       = variables_2d(i,j)%href
                sotypref     = variables_2d(i,j)%sotyp
                pressfcref   = variables_2d(i,j)%pressfc
                dswrfref     = variables_2d(i,j)%dswrf
                shtflref     = variables_2d(i,j)%shtfl
                tmpsfcref    = variables_2d(i,j)%tmpsfc
                tmp2mref     = variables_2d(i,j)%tmp2m
                spfh2mref    = variables_2d(i,j)%spfh2m
                hpblref      = variables_2d(i,j)%hpbl
                prate_averef = variables_2d(i,j)%prate_ave

! ... calculate wind speed from u and v
                ubzref   = sqrt((uref**2.0) + (vref**2.0))

! ... get scaled canopy model profile and sub-canopy layers
                zhc         = zk/hcmref
                cansublays  = floor(hcmref/modres)

! ... check for valid model vegetation types
                if (lu_opt .eq. 0 .or. lu_opt .eq. 1 ) then !VIIRS or MODIS
                    if (vtyperef .gt. 0 .and. vtyperef .le. 10 .or. vtyperef .eq. 12) then !VIIRS or MODIS types

! ... check for ssg_opt from user namelist
                        if (vtyperef .ge. 6 .and. vtyperef .le. 10) then !VIIRS/MODIS shrubs/savannas/grasses (SSG) type
                            if (ssg_opt .eq. 0) then !use GEDI inputs for SSG heights (not likely captured...)
                                hcmref = hcmref
                            else if (ssg_opt .eq. 1) then !user set constant shrubs/savannas/grasslands height
                                hcmref = ssg_set
                                !recalculate
                                zhc         = zk/hcmref
                                cansublays  = floor(hcmref/modres)
                            else
                                write(*,*)  'Wrong SSG_OPT choice of ', ssg_opt, &
                                    ' in namelist...exiting'
                                call exit(2)
                            end if
                        end if

! ... check for crop_opt from user namelist
                        if (vtyperef .eq. 12) then !VIIRS/MODIS crop type
                            if (crop_opt .eq. 0) then !use GEDI inputs for crop height (not likely captured...)
                                hcmref = hcmref
                            else if (crop_opt .eq. 1) then !user set constant crop height
                                hcmref = crop_set
                                !recalculate
                                zhc         = zk/hcmref
                                cansublays  = floor(hcmref/modres)
                            else
                                write(*,*)  'Wrong CROP_OPT choice of ', crop_opt, &
                                    ' in namelist...exiting'
                                call exit(2)
                            end if
                        end if

! ... check for contiguous canopy conditions at each model grid cell
                        if (hcmref .gt. fch_thresh .and. ffracref .gt. frt_thresh &
                            .and. lairef .gt. lai_thresh) then

! ... call canopy parameters to get canopy, fire info, and shape distribution parameters

                            call canopy_parm(vtyperef, hcmref, ffracref, lairef, &
                                pai_opt, pai_set, lu_opt, firetype, cdrag, &
                                pai, zcanmax, sigmau, sigma1)

! ... Choose between prescribed canopy/foliate shape profile or observed GEDI PAVD profile

                            if (pavd_opt .eq. 0) then
! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
                                call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                    fafraczInt)
!                                if (i .eq. 26 .and. j .eq. 26) then
!                                    print*,'prescribed shape function =', fafraczInt
!                                    print*,'zhc = ', zhc
!                                end if !test debug
                            else
! ... derive canopy/foliage distribution shape profile from interpolated GEDI PAVD profile - bottom up total in-canopy and fraction at z
!                                if (i .eq. 26 .and. j .eq. 26) then  !test debug
!                                    print*, 'Lat = ', variables_2d(i,j)%lat
!                                    print*, 'Lon = ', variables_2d(i,j)%lon
!                                    print*, 'VIIRS PAI = ', lairef+0.5
                                if (variables_2d(i,j)%lat .gt. (-1.0_rk*pavd_set) .and. &
                                    variables_2d(i,j)%lat .lt. pavd_set) then !use GEDI PAVD
                                    call canopy_pavd2fafrac(modres, hcmref, zhc, &
                                        variables_3d(i,j,:)%pavd, variables_1d%lev, fafraczInt)
!                                        print*, 'i = ', i, 'j = ', j
!                                        print*, 'fafraczInt(pavd) = ', fafraczInt
!                                        print*,'zhc = ', zhc
                                    !check if there is observed canopy height but no PAVD profile
                                    if (hcmref .gt. 0.0 .and. maxval(fafraczInt) .le. 0.0) then !revert to prescribed shape profile
                                        call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                            fafraczInt)
                                    end if
                                else !revert back to using prescribed shape profile
                                    call canopy_foliage(modlays, zhc, zcanmax, sigmau, sigma1, &
                                        fafraczInt)
                                end if
!                                end if !test debug
                            end if

! ... calculate zero-plane displacement height/hc and surface (soil+veg) roughness lengths/hc

                            call canopy_zpd(zhc(1:cansublays), fafraczInt(1:cansublays), &
                                ubzref, z0ghc, lambdars, cdrag, pai, hcmref, hgtref, &
                                z0ref, vtyperef, lu_opt, z0_opt, d_h, zo_h)

! ... calculate canopy radiation (sunlit and shade) profile

                            call canopy_fsun_clu( fafraczInt, lairef, cluref, cszref, fsun)

! ... calculate canopy leaf temperature (sun/shade) profile

                            call canopy_tleaf_lin(zk, hcmref, tmp2mref, fsun, &
                                tleaf_sun, tleaf_shade, tleaf_ave)

! ... calculate canopy Photosynthetic Photon Flux Density (PPFD) (sun/shade) profile

                            call canopy_ppfd_exp(zk, hcmref, dswrfref, lairef, fsun, &
                                ppfd_sun, ppfd_shade, ppfd_ave)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                            if (ifcanwind .or. ifcanwaf) then
                                if (rsl_opt .eq. 0) then
                                    do k=1, modlays
                                        call canopy_wind_most(hcmref, zk(k), fafraczInt(k), ubzref, &
                                            z0ghc, cdrag, pai, hgtref, d_h, zo_h, &
                                            lambdars, canBOT(k), canTOP(k), canWIND_3d(i,j,k))
                                    end do
                                else
                                    write(*,*) 'wrong namelist option = ', rsl_opt, 'only option = 0 right now'
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
                                                firetype, d_h, zo_h, canBOT(midflamepoint), &
                                                canTOP(midflamepoint), waf_2d(i,j))
                                        else
                                            write(*,*) 'warning...sub-canopy type fire, but flameh > FCH, setting WAF=1'
                                            waf_2d(i,j) = 1.0_rk
                                        end if
                                    end if
                                else  !grass/crops, above-canopy firetype
                                    if (flameh_2d(i,j) .gt. 0.0) then !flameh still must be > 0
                                        call canopy_waf(hcmref, lambdars, hgtref, flameh_2d(i,j), &
                                            firetype, d_h, zo_h, canBOT(midflamepoint), &
                                            canTOP(midflamepoint), waf_2d(i,j))
                                    end if
                                end if
                                if (waf_2d(i,j) .gt. 1.0_rk) then !Final check of WAF > 1, must be <=1
                                    waf_2d(i,j) = 1.0_rk
                                end if
                            end if

! ... user option to calculate in-canopy eddy diffusivities at height z
                            if (ifcaneddy) then
                                do k=1, modlays
                                    call canopy_eddyx(hcmref, zk(k), ustref, molref, Kz_3d(i,j,k))
                                end do
                            end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                            if (ifcanphot) then
                                if (cszref .ge. 0.0_rk) then !only calculate if cell isn't dark
                                    call canopy_phot(fafraczInt, &
                                        lairef, cluref, cszref, rjcf_3d(i,j,:))
                                end if
                            end if

!.......user option to calculate in-canopy leafage influence and assigning LAI as per timestep


                            if (leafage_opt .eq. 0) then
                                ! Initialize pastlai and currentlai based on current timestep
                                if (nn .eq. 1) then
                                    currentlai = lairef
                                    pastlai = currentlai
                                else
                                    pastlai = currentlai
                                    currentlai = lairef
                                end if

                                !!! Check if the lai_tstep is greater than time_intvl
                                if (lai_tstep .ge. time_intvl) then
                                    tsteplai = lai_tstep/86400.0_rk
                                else
                                    WRITE (*, *) "Error: Input LAI time step cannot be less than model time step...exiting!!!"
                                    CALL EXIT(1)
                                endif
                            end if !leafage_opt = 0 end


! ... user option to calculate in-canopy biogenic emissions
                            if (ifcanbio) then
                                if (cszref .ge. 0.0_rk .and. dswrfref .gt. 0.0_rk &
                                    .and. cluref .gt. 0.0_rk) then
                                    !ISOP
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        1, emi_isop_3d(i,j,:))
                                    !MYRC
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        2, emi_myrc_3d(i,j,:))
                                    !SABI
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        3, emi_sabi_3d(i,j,:))
                                    !LIMO
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        4, emi_limo_3d(i,j,:))
                                    !CARE
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        5, emi_care_3d(i,j,:))
                                    !OCIM
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        6, emi_ocim_3d(i,j,:))
                                    !BPIN
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        7, emi_bpin_3d(i,j,:))
                                    !APIN
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        8, emi_apin_3d(i,j,:))
                                    !MONO
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        9, emi_mono_3d(i,j,:))
                                    !FARN
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        10, emi_farn_3d(i,j,:))
                                    !CARY
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        11, emi_cary_3d(i,j,:))
                                    !SESQ
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        12, emi_sesq_3d(i,j,:))
                                    !MBOL
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        13, emi_mbol_3d(i,j,:))
                                    !METH
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        14, emi_meth_3d(i,j,:))
                                    !ACET
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        15, emi_acet_3d(i,j,:))
                                    !CO
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        16, emi_co_3d(i,j,:))
                                    !BIDI VOC
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        17, emi_bvoc_3d(i,j,:))
                                    !Stress VOC
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        18, emi_svoc_3d(i,j,:))
                                    !Other VOC
                                    call canopy_bio(zk, fafraczInt, hcmref, &
                                        lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                        tleaf_ave, tmp2mref, &
                                        lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                        leafage_opt, pastlai, currentlai, tsteplai,  &
                                        19, emi_ovoc_3d(i,j,:))
                                end if
                            end if

                        end if !Contiguous Canopy

                    else
!                        write(*,*)  'Warning VIIRS/MODIS VTYPE ', vtyperef, ' is not supported...continue'
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
            call canopy_calcdx(dx_opt, dx_set, nlat, nlon, variables%lat, &
                variables%lon, dx)
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
            hcmref       = variables(loc)%fh
            uref         = variables(loc)%ugrd10m
            vref         = variables(loc)%vgrd10m
            cluref       = variables(loc)%clu
            lairef       = variables(loc)%lai
            vtyperef     = variables(loc)%vtype
            ffracref     = variables(loc)%ffrac
            ustref       = variables(loc)%fricv
            cszref       = variables(loc)%csz
            z0ref        = variables(loc)%sfcr
            molref       = variables(loc)%mol
            frpref       = variables(loc)%frp
            hgtref       = variables(loc)%href
            sotypref     = variables(loc)%sotyp
            pressfcref   = variables(loc)%pressfc
            dswrfref     = variables(loc)%dswrf
            shtflref     = variables(loc)%shtfl
            tmpsfcref    = variables(loc)%tmpsfc
            tmp2mref     = variables(loc)%tmp2m
            spfh2mref    = variables(loc)%spfh2m
            hpblref      = variables(loc)%hpbl
            prate_averef = variables(loc)%prate_ave

! ... calculate wind speed from u and v
            ubzref   = sqrt((uref**2.0) + (vref**2.0))

! ... get scaled canopy model profile and sub-canopy layers
            zhc         = zk/hcmref
            cansublays  = floor(hcmref/modres)

! ... check for valid model vegetation types
            if (lu_opt .eq. 0 .or. lu_opt .eq. 1 ) then !VIIRS or MODIS
                if (vtyperef .gt. 0 .and. vtyperef .le. 10 .or. vtyperef .eq. 12) then !VIIRS or MODIS types

! ... check for ssg_opt from user namelist
                    if (vtyperef .ge. 6 .and. vtyperef .le. 10) then !VIIRS/MODIS shrubs/savannas/grasses (SSG) type
                        if (ssg_opt .eq. 0) then !use GEDI inputs for SSG heights (not likely captured...)
                            hcmref = hcmref
                        else if (ssg_opt .eq. 1) then !user set constant shrubs/savannas/grasslands height
                            hcmref = ssg_set
                            !recalculate
                            zhc         = zk/hcmref
                            cansublays  = floor(hcmref/modres)
                        else
                            write(*,*)  'Wrong SSG_OPT choice of ', ssg_opt, &
                                ' in namelist...exiting'
                            call exit(2)
                        end if
                    end if

! ... check for crop_opt from user namelist
                    if (vtyperef .eq. 12) then !VIIRS/MODIS crop type
                        if (crop_opt .eq. 0) then !use GEDI inputs for crop height
                            hcmref = hcmref
                        else if (crop_opt .eq. 1) then !user set constant crop height
                            hcmref = crop_set
                            !recalculate
                            zhc         = zk/hcmref
                            cansublays  = floor(hcmref/modres)
                        else
                            write(*,*)  'Wrong CROP_OPT choice of ', crop_opt, &
                                ' in namelist...exiting'
                            call exit(2)
                        end if
                    end if

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
                            ubzref, z0ghc, lambdars, cdrag, pai, hcmref, hgtref, &
                            z0ref, vtyperef, lu_opt, z0_opt, d_h, zo_h)

! ... calculate canopy radiation (sunlit and shade) profile

                        call canopy_fsun_clu( fafraczInt, lairef, cluref, cszref, fsun)

! ... calculate canopy leaf temperature (sun/shade) profile

                        call canopy_tleaf_lin(zk, hcmref, tmp2mref, fsun, &
                            tleaf_sun, tleaf_shade, tleaf_ave)

! ... calculate canopy Photosynthetic Photon Flux Density (PPFD) (sun/shade) profile

                        call canopy_ppfd_exp(zk, hcmref, dswrfref, lairef, fsun, &
                            ppfd_sun, ppfd_shade, ppfd_ave)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                        if (ifcanwind .or. ifcanwaf) then
                            if (rsl_opt .eq. 0) then
                                do k=1, modlays
                                    call canopy_wind_most(hcmref, zk(k), fafraczInt(k), ubzref, &
                                        z0ghc, cdrag, pai, hgtref, d_h, zo_h, &
                                        lambdars, canBOT(k), canTOP(k), canWIND(loc,k))
                                end do
                            else
                                write(*,*) 'wrong namelist option = ', rsl_opt, 'only option = 0 right now'
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
                                            firetype, d_h, zo_h, canBOT(midflamepoint), &
                                            canTOP(midflamepoint), waf(loc))
                                    else
                                        write(*,*) 'warning...sub-canopy type fire, but flameh > FCH, setting WAF=1'
                                        waf(loc) = 1.0_rk
                                    end if
                                end if
                            else  !grass/crops, above-canopy firetype
                                if (flameh(loc) .gt. 0.0) then !flameh still must be > 0
                                    call canopy_waf(hcmref, lambdars, hgtref, flameh(loc), &
                                        firetype, d_h, zo_h, canBOT(midflamepoint), &
                                        canTOP(midflamepoint), waf(loc))
                                end if
                            end if
                            if (waf(loc) .gt. 1.0_rk) then !Final check of WAF > 1, must be <=1
                                waf(loc) = 1.0_rk
                            end if
                        end if

! ... user option to calculate in-canopy eddy diffusivities at height z
                        if (ifcaneddy) then
                            do k=1, modlays
                                call canopy_eddyx(hcmref, zk(k), ustref, molref, Kz(loc, k))
                            end do
                        end if

! ... user option to calculate in-canopy eddy photolysis attenuation at height z
                        if (ifcanphot) then
                            if (cszref .ge. 0.0_rk) then !only calculate if cell isn't dark
                                call canopy_phot(fafraczInt, &
                                    lairef, cluref, cszref, rjcf(loc, :))
                            end if
                        end if

! ... user option to calculate in-canopy biogenic emissions
                        if (ifcanbio) then
                            if (cszref .ge. 0.0_rk .and. dswrfref .gt. 0.0_rk &
                                .and. cluref .gt. 0.0_rk) then
                                !ISOP
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    1, emi_isop(loc,:))
                                !MYRC
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    2, emi_myrc(loc,:))
                                !SABI
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    3, emi_sabi(loc,:))
                                !LIMO
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    4, emi_limo(loc,:))
                                !CARE
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    5, emi_care(loc,:))
                                !OCIM
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    6, emi_ocim(loc,:))
                                !BPIN
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    7, emi_bpin(loc,:))
                                !APIN
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    8, emi_apin(loc,:))
                                !MONO
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    9, emi_mono(loc,:))
                                !FARN
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    10, emi_farn(loc,:))
                                !CARY
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    11, emi_cary(loc,:))
                                !SESQ
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    12, emi_sesq(loc,:))
                                !MBOL
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    13, emi_mbol(loc,:))
                                !METH
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    14, emi_meth(loc,:))
                                !ACET
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    15, emi_acet(loc,:))
                                !CO
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    16, emi_co(loc,:))
                                !BIDI VOC
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    17, emi_bvoc(loc,:))
                                !Stress VOC
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    18, emi_svoc(loc,:))
                                !Other VOC
                                call canopy_bio(zk, fafraczInt, hcmref, &
                                    lairef, fsun, ppfd_sun, ppfd_shade, tleaf_sun, tleaf_shade, &
                                    tleaf_ave, tmp2mref, &
                                    lu_opt, vtyperef, modres, bio_cce, biovert_opt, co2_opt, co2_set, &
                                    leafage_opt, pastlai, currentlai, tsteplai,  &
                                    19, emi_ovoc(loc,:))
                            end if
                        end if

                    end if !Contiguous Canopy

                else
!                    write(*,*)  'Warning VIIRS/MODIS VTYPE ', vtyperef, ' is not supported...continue'
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
