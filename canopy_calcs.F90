
SUBROUTINE canopy_calcs

!-------------------------------------------------------------------------------
! Name:     Contains main canopy calculations dependent on canopy conditions
! Purpose:  Contains main canopy calculations dependent on canopy conditions
! Revised:  06 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------


    use canopy_coord_mod   !main canopy coordinate descriptions
    use canopy_canopts_mod !main canopy option descriptions
    use canopy_canmet_mod  !main canopy met/sfc input descriptions
    use canopy_canvars_mod !main canopy variables descriptions

    use canopy_driver_mod  !main canopy components/submodules

    IMPLICIT NONE

    !Local variables
    integer i,loc


    write(*,*)  'Calculating Canopy Parameters'
    write(*,*)  '-------------------------------'

    if (ifcanwind) then !only calculate if canopy wind option
        call canopy_calcdx(dx_opt, dx_set, nlat, nlon, variables%lat, &
            variables%lon, dx)
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
            end if   !Vegetation types
        else
            write(*,*)  'Wrong LU_OPT choice of ', lu_opt, ' in namelist...exiting'
            call exit(2)

        end if       !Landuse Options
    end do


END SUBROUTINE canopy_calcs
