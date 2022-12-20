
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

    IMPLICIT NONE

    !Local variables
    integer i,loc
    real(rk) hgtref !local reference height from namelist or file array

    write(*,*)  'Calculating Canopy Parameters'
    write(*,*)  '-------------------------------'

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
        hgtref   = variables(loc)%href

! ... get scaled canopy model profile and sub-canopy layers
        zhc         = zk/hcmref
        cansublays  = floor(hcmref/modres)

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
                        ubzref, z0ghc, lamdars, rsl_opt, cdrag, pai, hcmref, hgtref, &
                        z0ref, vtyperef, lu_opt, z0_opt, d_h, zo_h)

! ... user option to calculate in-canopy wind speeds at height z and midflame WAF

                    if (ifcanwind .or. ifcanwaf) then
                        do i=1, modlays
                            call canopy_wind(hcmref, zk(i), fafraczInt(i), ubzref, &
                                z0ghc, cdrag, pai, hgtref, d_h, zo_h, molref, &
                                rsl_opt, canBOT(i), canTOP(i), canWIND(i, loc))
                        end do
! ... determine midflamepoint and flame height from user or FRP calculation
                        call canopy_flameh(flameh_opt, flameh_set, dx(loc), modres, &
                            frpref, midflamepoint, flameh)

                        if (flameh .gt. 0.0) then !only calculate WAF when flameh > 0
                            call canopy_waf(hcmref, lamdars, rsl_opt, hgtref, flameh, &
                                firetype, d_h, zo_h, canBOT(midflamepoint), &
                                canTOP(midflamepoint), waf(loc))
                        end if
                    end if

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

            else
                write(*,*)  'Warning VIIRS VTYPE ', vtyperef, ' is not supported...continue'
            end if   !Vegetation types
        else
            write(*,*)  'Wrong LU_OPT choice of ', lu_opt, ' in namelist...exiting'
            call exit(2)

        end if       !Landuse Options
    end do


END SUBROUTINE canopy_calcs
