program canopy_driver
!
!  The main driver to run the canopy applications
!
!  History:
!    Prototype: Patrick C. Campbell, 06/2022
!
!-------------------------------------------------------------
    use canopy_coord_mod   !main canopy coordinate descriptions
    use canopy_canopts_mod !main canopy option descriptions
    use canopy_canmet_mod  !main canopy met/sfc input descriptions
    use canopy_canvars_mod !main canopy variables descriptions
    use canopy_files_mod   !main canopy input files

    use canopy_calcdx_mod  !main grid dx calculation
    use canopy_parm_mod    !main canopy parameters
    use canopy_foliage_mod !main canopy foliage distribution
    use canopy_zpd_mod     !main displacement height model
    use canopy_wind_mod    !main canopy wind model
    use canopy_flameh_mod  !main flame height model
    use canopy_waf_mod     !main Wind Adjustment Factor (WAF) model
    use canopy_eddyx_mod   !main canopy eddy diffusivities
    use canopy_phot_mod    !main canopy photolysis attenuation

    implicit none

!Local variables
    integer i,loc

!-------------------------------------------------------------------------------
! Read user options from namelist.
!-------------------------------------------------------------------------------

    call canopy_readnml

!-------------------------------------------------------------------------------
! Allocate necessary variables.
!-------------------------------------------------------------------------------

    call canopy_alloc

!-------------------------------------------------------------------------------
! Read met/sfc gridded model input file.
!-------------------------------------------------------------------------------

! ... TODO:  NL condition for data read from txt or netcdf 1D or 2D

    call canopy_read_txt(file_vars(1))

! ... TODO:  add option to read met/sfc input variables from 1D ncf file

    !call canopy_read_ncf_1D(file_vars(1))

! ... TODO:  add option to read met/sfc input variables from 2D ncf file
    !call canopy_read_ncf_2D(file_vars(1))

!-------------------------------------------------------------------------------
! Initialize canopy grid cell dependent variables
!-------------------------------------------------------------------------------

! ... TODO:  add canopy_grid initialize subroutine

! ... initialize ***grid cell_only*** dependent variables
    dx   = 0.0
    waf  = 1.0

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
        fafraczInt        = 0.0
        canTOP            = 1.0
        canBOT            = 1.0

! ... initialize ***grid cell_and_canopy*** profile dependent variables
        canWIND(:,loc)    = ubzref    !initialize to above canopy wind
        Kz(:,loc)         = -999.0    !initialize to missing
        rjcf(:,loc)       = 1.0       !initialize to above canopy rjcf=1

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

    call canopy_dealloc

end program canopy_driver
