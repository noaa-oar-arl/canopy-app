
SUBROUTINE canopy_calcs

!-------------------------------------------------------------------------------
! Name:     Contains main canopy calculations dependent on canopy conditions
! Purpose:  Contains main canopy calculations dependent on canopy conditions
! Revised:  06 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    use canopy_const_mod, ONLY: rk      !constants for canopy models
    use canopy_coord_mod   !main canopy coordinate descriptions
    use canopy_canopts_mod !main canopy option descriptions
    use canopy_canmet_mod  !main canopy met/sfc input descriptions
    use canopy_canvars_mod !main canopy variables descriptions
    use canopy_utils_mod   !main canopy utilities
    use canopy_dxcalc_mod  !main canopy dx calculation
    use canopy_profile_mod !main canopy foliage profile routines
    use canopy_wind_mod    !main canopy components/submodules
    use canopy_waf_mod
    use canopy_phot_mod
    use canopy_eddy_mod
    use canopy_bioemi_mod

    IMPLICIT NONE

    !Local variables
    integer i,j,k,loc

    write(*,*)  'Calculating Canopy Parameters'
    write(*,*)  '-------------------------------'

    if (infmt_opt .eq. 0) then !Input format is 2D then output is then 2D

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

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                            if (ifcanwind .or. ifcanwaf) then
                                do k=1, modlays
                                    call canopy_wind(hcmref, zk(k), fafraczInt(k), ubzref, &
                                        z0ghc, cdrag, pai, hgtref, d_h, zo_h, molref, &
                                        rsl_opt, lambdars, canBOT(k), canTOP(k), canWIND_3d(i,j,k))
                                end do

! ... determine midflamepoint and flame height from user or FRP calculation
                                call canopy_flameh(flameh_opt, flameh_set, dx_2d(i,j), modres, &
                                    frpref, midflamepoint, flameh_2d(i,j))

                                if (flameh_2d(i,j) .gt. 0.0 .and. flameh_2d(i,j) .le. hcmref) then
                                    !only calculate WAF when flameh > 0 and <= FH
                                    call canopy_waf(hcmref, lambdars, hgtref, flameh_2d(i,j), &
                                        firetype, d_h, zo_h, canBOT(midflamepoint), &
                                        canTOP(midflamepoint), waf_2d(i,j))
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

! ... user option to calculate in-canopy biogenic emissions
                            if (ifcanbio) then
                               if (cszref .ge. 0.0_rk .and. dswrfref .gt. 0.0_rk) then 
                                          !ISOP
                                          call canopy_bio(zk, fafraczInt, hcmref, &
                                          lairef, cluref, cszref, dswrfref, tmp2mref, &
                                          lu_opt, vtyperef, 0, emi_isop_3d(i,j,:))
                                          !MYRC
                                          call canopy_bio(zk, fafraczInt, hcmref, &
                                          lairef, cluref, cszref, dswrfref, tmp2mref, &
                                          lu_opt, vtyperef, 1, emi_myrc_3d(i,j,:))
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

    else if (infmt_opt .eq. 1) then !Input format is 1D and output must be 1D

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

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                        if (ifcanwind .or. ifcanwaf) then
                            do k=1, modlays
                                call canopy_wind(hcmref, zk(k), fafraczInt(k), ubzref, &
                                    z0ghc, cdrag, pai, hgtref, d_h, zo_h, molref, &
                                    rsl_opt, lambdars, canBOT(k), canTOP(k), canWIND(loc, k))
                            end do
! ... determine midflamepoint and flame height from user or FRP calculation
                            call canopy_flameh(flameh_opt, flameh_set, dx(loc), modres, &
                                frpref, midflamepoint, flameh(loc))

                            if (flameh(loc) .gt. 0.0 .and. flameh(loc) .le. hcmref) then
                                !only calculate WAF when flameh > 0
                                call canopy_waf(hcmref, lambdars, hgtref, flameh(loc), &
                                    firetype, d_h, zo_h, canBOT(midflamepoint), &
                                    canTOP(midflamepoint), waf(loc))
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
                               if (cszref .ge. 0.0_rk .and. dswrfref .gt. 0.0_rk) then 
                                          !ISOP
                                          call canopy_bio(zk, fafraczInt, hcmref, &
                                          lairef, cluref, cszref, dswrfref, tmp2mref, &
                                          lu_opt, vtyperef, 0, emi_isop(loc,:))
                                          !MYRC
                                          call canopy_bio(zk, fafraczInt, hcmref, &
                                          lairef, cluref, cszref, dswrfref, tmp2mref, &
                                          lu_opt, vtyperef, 1, emi_myrc(loc,:)) 
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

    end if !Input Format (1D or 2D)

END SUBROUTINE canopy_calcs
