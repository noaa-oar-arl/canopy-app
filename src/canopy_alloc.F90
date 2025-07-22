
!> \file canopy_alloc.F90
!! \brief Memory allocation subroutine for canopy model arrays
!! \author P.C. Campbell
!! \date 03 Oct 2022
!! \version 1.0

!> \defgroup MemoryManagement Memory Management
!! \brief Routines    !> \defgroup Daily3DAverages Daily Average 3D Arrays
!! \brief 3D versions of daily average arrays
!! \ingroup MemoryManagement
!! \{
!> \brief 24-hour average sunlit leaf temperature 3D (K)locating and deallocating memory for canopy model arrays
!! \{

!> \brief Allocate all necessary arrays for canopy model inputs and outputs
!!
!! This subroutine allocates memory for all arrays used in the canopy model
!! including input variables, canopy distribution variables, wind fields,
!! diffusivity profiles, photolysis correction factors, biogenic emissions,
!! and gas dry deposition arrays. Memory allocation is conditional based on
!! user-selected options.
!!
!! \author P.C. Campbell
!! \date 03 Oct 2022
!! \version 1.0
!!
!! \par Revision History:
!! - 03 Oct 2022: Original version (P.C. Campbell)
!!
!! \par References:
!! - Canopy-App Documentation
!!
!! \ingroup MemoryManagement
SUBROUTINE canopy_alloc

    !> \brief Module for canopy model options and parameters
    USE canopy_canopts_mod
    !> \brief Module for coordinate and domain information
    USE canopy_coord_mod
    !> \brief Module for canopy meteorological input variables
    USE canopy_canmet_mod
    !> \brief Module for canopy model variables and arrays
    USE canopy_canvars_mod

    IMPLICIT NONE

    !> \defgroup InputVariables Input Variable Arrays
    !! \brief Arrays for storing meteorological and surface input data
    !! \ingroup MemoryManagement
    !! \{

    !> \brief 1D array of input variables for each grid cell
    if(.not.allocated(variables))     allocate(variables(nlat*nlon))
    !> \brief 2D array of input variables organized by longitude and latitude
    if(.not.allocated(variables_2d))  allocate(variables_2d(nlon,nlat))

    !> \brief Conditional allocation for 3D variable arrays
    if (var3d_opt .eq. 1) then
        !> \brief 3D array for variables with vertical dimension
        if(.not.allocated(variables_3d))  allocate(variables_3d(nlon,nlat,var3d_set))
        !> \brief 1D array for vertical level information
        if(.not.allocated(variables_1d))  allocate(variables_1d(var3d_set))
        !> \brief 1D array for canopy-specific variables
        if(.not.allocated(variables_can)) allocate(variables_can(nlat*nlon))
        !> \brief Plant Area Volume Density reference profile
        if(.not.allocated(pavdref))       allocate(pavdref(var3d_set))
        !> \brief Reference vertical levels array
        if(.not.allocated(levref))        allocate(levref(var3d_set))
        !> \brief Working array for PAVD values
        if(.not.allocated(pavd_arr))      allocate(pavd_arr(var3d_set))
        !> \brief Working array for level values
        if(.not.allocated(lev_arr))       allocate(lev_arr(var3d_set))
    end if
    !> \}

    !> \defgroup CanopyDistribution Canopy Distribution Arrays
    !! \brief Arrays for storing canopy structure and distribution information
    !! \ingroup MemoryManagement
    !! \{

    !> \brief Normalized height within canopy (z/hc)
    if(.not.allocated(zhc))                allocate(zhc(modlays))
    !> \brief Fractional cumulative leaf area index profile
    if(.not.allocated(fafraczInt))         allocate(fafraczInt(modlays))
    !> \brief Sunlit leaf fraction profile
    if(.not.allocated(fsun))               allocate(fsun(modlays))
    !> \brief Sunlit leaf temperature profile (K)
    if(.not.allocated(tleaf_sun))          allocate(tleaf_sun(modlays))
    !> \brief Shaded leaf temperature profile (K)
    if(.not.allocated(tleaf_shade))        allocate(tleaf_shade(modlays))
    !> \brief Average leaf temperature profile (K)
    if(.not.allocated(tleaf_ave))          allocate(tleaf_ave(modlays))
    !> \brief Sunlit leaf photosynthetic photon flux density (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_sun))           allocate(ppfd_sun(modlays))
    !> \brief Shaded leaf photosynthetic photon flux density (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_shade))         allocate(ppfd_shade(modlays))
    !> \brief Average photosynthetic photon flux density (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_ave))           allocate(ppfd_ave(modlays))
    !> \brief Leaf area density profile 1D (m² m⁻³)
    if(.not.allocated(lad))                allocate(lad(nlat*nlon,modlays))
    !> \brief Leaf area density profile 3D (m² m⁻³)
    if(.not.allocated(lad_3d))             allocate(lad_3d(nlon,nlat,modlays))
    !> \brief Surface roughness length normalized by canopy height
    if(.not.allocated(zo_h))               allocate(zo_h(nlat*nlon))
    !> \brief Surface roughness length normalized by canopy height (2D)
    if(.not.allocated(zo_h_2d))            allocate(zo_h_2d(nlon,nlat))
    !> \brief Displacement height normalized by canopy height
    if(.not.allocated(d_h))                allocate(d_h(nlat*nlon))
    !> \brief Displacement height normalized by canopy height (2D)
    if(.not.allocated(d_h_2d))             allocate(d_h_2d(nlon,nlat))
    !> \brief Air temperature profile within canopy (K)
    if(.not.allocated(tka))                allocate(tka(nlat*nlon,modlays))
    !> \brief Air temperature profile within canopy 3D (K)
    if(.not.allocated(tka_3d))             allocate(tka_3d(nlon,nlat,modlays))
    !> \brief Air pressure profile within canopy (Pa)
    if(.not.allocated(pressa))             allocate(pressa(nlat*nlon,modlays))
    !> \brief Air pressure profile within canopy 3D (Pa)
    if(.not.allocated(pressa_3d))          allocate(pressa_3d(nlon,nlat,modlays))
    !> \brief Relative humidity profile within canopy (%)
    if(.not.allocated(relhuma))            allocate(relhuma(nlat*nlon,modlays))
    !> \brief Relative humidity profile within canopy 3D (%)
    if(.not.allocated(relhuma_3d))         allocate(relhuma_3d(nlon,nlat,modlays))
    !> \brief Specific humidity profile within canopy (kg kg⁻¹)
    if(.not.allocated(spechuma))           allocate(spechuma(nlat*nlon,modlays))
    !> \brief Specific humidity profile within canopy 3D (kg kg⁻¹)
    if(.not.allocated(spechuma_3d))        allocate(spechuma_3d(nlon,nlat,modlays))

    !> \defgroup TemporalArrays Temporal History Arrays
    !! \brief Arrays for storing temporal histories when hist_opt=1
    !! \ingroup MemoryManagement
    !! \{
    if (hist_opt .eq. 1) then
        !> \brief 24-hour sunlit leaf temperature history (25 hours, K)
        if(.not.allocated(tleaf_sun24_tmp))    allocate(tleaf_sun24_tmp(25,nlat*nlon,modlays))
        !> \brief 24-hour shaded leaf temperature history (25 hours, K)
        if(.not.allocated(tleaf_shade24_tmp))  allocate(tleaf_shade24_tmp(25,nlat*nlon,modlays))
        !> \brief 24-hour average leaf temperature history (25 hours, K)
        if(.not.allocated(tleaf_ave24_tmp))    allocate(tleaf_ave24_tmp(25,nlat*nlon,modlays))
        !> \brief 24-hour sunlit PPFD history (25 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_sun24_tmp))     allocate(ppfd_sun24_tmp(25,nlat*nlon,modlays))
        !> \brief 24-hour shaded PPFD history (25 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_shade24_tmp))   allocate(ppfd_shade24_tmp(25,nlat*nlon,modlays))
        !> \brief 240-hour sunlit leaf temperature history (241 hours, K)
        if(.not.allocated(tleaf_sun240_tmp))   allocate(tleaf_sun240_tmp(241,nlat*nlon,modlays))
        !> \brief 240-hour shaded leaf temperature history (241 hours, K)
        if(.not.allocated(tleaf_shade240_tmp)) allocate(tleaf_shade240_tmp(241,nlat*nlon,modlays))
        !> \brief 240-hour average leaf temperature history (241 hours, K)
        if(.not.allocated(tleaf_ave240_tmp))   allocate(tleaf_ave240_tmp(241,nlat*nlon,modlays))
        !> \brief 240-hour sunlit PPFD history (241 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_sun240_tmp))    allocate(ppfd_sun240_tmp(241,nlat*nlon,modlays))
        !> \brief 240-hour shaded PPFD history (241 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_shade240_tmp))  allocate(ppfd_shade240_tmp(241,nlat*nlon,modlays))
        !> \brief 24-hour 2m temperature reference history (25 hours, K)
        if(.not.allocated(tmp2mref_tmp))       allocate(tmp2mref_tmp(25,nlat*nlon))
        !> \brief 24-hour reference wind speed history (25 hours, m s⁻¹)
        if(.not.allocated(ubzref_tmp))         allocate(ubzref_tmp(25,nlat*nlon))
    end if

    !> \defgroup DailyAverages Daily Average Arrays
    !! \brief Arrays for storing daily average values
    !! \ingroup MemoryManagement
    !! \{
    !> \brief 24-hour average sunlit leaf temperature (K)
    if(.not.allocated(tleaf_sun24))        allocate(tleaf_sun24(nlat*nlon,modlays))
    !> \brief 24-hour average shaded leaf temperature (K)
    if(.not.allocated(tleaf_shade24))      allocate(tleaf_shade24(nlat*nlon,modlays))
    !> \brief 24-hour average leaf temperature (K)
    if(.not.allocated(tleaf_ave24))        allocate(tleaf_ave24(nlat*nlon,modlays))
    !> \brief 24-hour average sunlit PPFD (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_sun24))         allocate(ppfd_sun24(nlat*nlon,modlays))
    !> \brief 24-hour average shaded PPFD (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_shade24))       allocate(ppfd_shade24(nlat*nlon,modlays))
    !> \brief 240-hour average sunlit leaf temperature (K)
    if(.not.allocated(tleaf_sun240))       allocate(tleaf_sun240(nlat*nlon,modlays))
    !> \brief 240-hour average shaded leaf temperature (K)
    if(.not.allocated(tleaf_shade240))     allocate(tleaf_shade240(nlat*nlon,modlays))
    !> \brief 240-hour average leaf temperature (K)
    if(.not.allocated(tleaf_ave240))       allocate(tleaf_ave240(nlat*nlon,modlays))
    !> \brief 240-hour average sunlit PPFD (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_sun240))        allocate(ppfd_sun240(nlat*nlon,modlays))
    !> \brief 240-hour average shaded PPFD (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_shade240))      allocate(ppfd_shade240(nlat*nlon,modlays))

    !> \brief Daily maximum 2m temperature (K)
    if(.not.allocated(daily_maxt2m))       allocate(daily_maxt2m(nlat*nlon))
    !> \brief Daily minimum 2m temperature (K)
    if(.not.allocated(daily_mint2m))       allocate(daily_mint2m(nlat*nlon))
    !> \brief Daily maximum 10m wind speed (cm s-1)
    if(.not.allocated(daily_maxws10m))     allocate(daily_maxws10m(nlat*nlon))

    !> \defgroup Temporal3DArrays Temporal History 3D Arrays
    !! \brief 3D versions of temporal history arrays when hist_opt=1
    !! \ingroup MemoryManagement
    !! \{
    if (hist_opt .eq. 1) then
        !> \brief 24-hour sunlit leaf temperature 3D history (25 hours, K)
        if(.not.allocated(tleaf_sun24_tmp_3d))    allocate(tleaf_sun24_tmp_3d(25,nlon,nlat,modlays))
        !> \brief 24-hour shaded leaf temperature 3D history (25 hours, K)
        if(.not.allocated(tleaf_shade24_tmp_3d))  allocate(tleaf_shade24_tmp_3d(25,nlon,nlat,modlays))
        !> \brief 24-hour average leaf temperature 3D history (25 hours, K)
        if(.not.allocated(tleaf_ave24_tmp_3d))    allocate(tleaf_ave24_tmp_3d(25,nlon,nlat,modlays))
        !> \brief 24-hour sunlit PPFD 3D history (25 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_sun24_tmp_3d))     allocate(ppfd_sun24_tmp_3d(25,nlon,nlat,modlays))
        !> \brief 24-hour shaded PPFD 3D history (25 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_shade24_tmp_3d))   allocate(ppfd_shade24_tmp_3d(25,nlon,nlat,modlays))
        !> \brief 240-hour sunlit leaf temperature 3D history (241 hours, K)
        if(.not.allocated(tleaf_sun240_tmp_3d))   allocate(tleaf_sun240_tmp_3d(241,nlon,nlat,modlays))
        !> \brief 240-hour shaded leaf temperature 3D history (241 hours, K)
        if(.not.allocated(tleaf_shade240_tmp_3d)) allocate(tleaf_shade240_tmp_3d(241,nlon,nlat,modlays))
        !> \brief 240-hour average leaf temperature 3D history (241 hours, K)
        if(.not.allocated(tleaf_ave240_tmp_3d))   allocate(tleaf_ave240_tmp_3d(241,nlon,nlat,modlays))
        !> \brief 240-hour sunlit PPFD 3D history (241 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_sun240_tmp_3d))    allocate(ppfd_sun240_tmp_3d(241,nlon,nlat,modlays))
        !> \brief 240-hour shaded PPFD 3D history (241 hours, µmol m⁻² s⁻¹)
        if(.not.allocated(ppfd_shade240_tmp_3d))  allocate(ppfd_shade240_tmp_3d(241,nlon,nlat,modlays))
        !> \brief 24-hour 2m temperature reference 3D history (25 hours, K)
        if(.not.allocated(tmp2mref_tmp_3d))       allocate(tmp2mref_tmp_3d(25,nlon,nlat))
        !> \brief 24-hour reference wind speed 3D history (25 hours, m s⁻¹)
        if(.not.allocated(ubzref_tmp_3d))         allocate(ubzref_tmp_3d(25,nlon,nlat))
    end if

    !>
    !! \defgroup Daily3DAverages Daily Average 3D Arrays
    !! \brief 3D versions of daily average arrays
    !! \ingroup MemoryManagement
    !! \{
    !!
    !<
    !> \brief 24-hour average sunlit leaf temperature 3D (K)
    if(.not.allocated(tleaf_sun24_3d))        allocate(tleaf_sun24_3d(nlon,nlat,modlays))
    !> \brief 24-hour average shaded leaf temperature 3D (K)
    if(.not.allocated(tleaf_shade24_3d))      allocate(tleaf_shade24_3d(nlon,nlat,modlays))
    !> \brief 24-hour average leaf temperature 3D (K)
    if(.not.allocated(tleaf_ave24_3d))        allocate(tleaf_ave24_3d(nlon,nlat,modlays))
    !> \brief 24-hour average sunlit PPFD 3D (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_sun24_3d))         allocate(ppfd_sun24_3d(nlon,nlat,modlays))
    !> \brief 24-hour average shaded PPFD 3D (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_shade24_3d))       allocate(ppfd_shade24_3d(nlon,nlat,modlays))
    !> \brief 240-hour average sunlit leaf temperature 3D (K)
    if(.not.allocated(tleaf_sun240_3d))       allocate(tleaf_sun240_3d(nlon,nlat,modlays))
    !> \brief 240-hour average shaded leaf temperature 3D (K)
    if(.not.allocated(tleaf_shade240_3d))     allocate(tleaf_shade240_3d(nlon,nlat,modlays))
    !> \brief 240-hour average leaf temperature 3D (K)
    if(.not.allocated(tleaf_ave240_3d))       allocate(tleaf_ave240_3d(nlon,nlat,modlays))
    !> \brief 240-hour average sunlit PPFD 3D (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_sun240_3d))        allocate(ppfd_sun240_3d(nlon,nlat,modlays))
    !> \brief 240-hour average shaded PPFD 3D (µmol m⁻² s⁻¹)
    if(.not.allocated(ppfd_shade240_3d))      allocate(ppfd_shade240_3d(nlon,nlat,modlays))

    !> \brief Daily maximum 2m temperature 3D (K)
    if(.not.allocated(daily_maxt2m_2d))       allocate(daily_maxt2m_2d(nlon,nlat))
    !> \brief Daily minimum 2m temperature 3D (K)
    if(.not.allocated(daily_mint2m_2d))       allocate(daily_mint2m_2d(nlon,nlat))
    !> \brief Daily maximum 10m wind speed 3D (cm s-1)
    if(.not.allocated(daily_maxws10m_2d))     allocate(daily_maxws10m_2d(nlon,nlat))
    !> \}

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Wind Outputs
!-------------------------------------------------------------------------------

    !> \defgroup CanopyWind Canopy Wind Arrays
    !! \brief Arrays for canopy wind speed profiles and wildfire applications
    !! \ingroup MemoryManagement
    !! \{
    !!
    !<
    if (ifcanwind .or. ifcanwaf) then
        write(*,*)  'Canopy wind and/or WAF option selected'
        write(*,*)  '-------------------------------'
        !> \brief Bottom height of canopy layers (m)
        if(.not.allocated(canBOT))        allocate(canBOT(modlays))
        !> \brief Top height of canopy layers (m)
        if(.not.allocated(canTOP))        allocate(canTOP(modlays))
        !> \brief Wind speed within canopy layers (m s-1)
        if(.not.allocated(canWIND))       allocate(canWIND(nlat*nlon,modlays))
        !> \brief Wind speed within canopy layers 3D (m s-1)
        if(.not.allocated(canWIND_3d))    allocate(canWIND_3d(nlon,nlat,modlays))
        !> \brief Grid cell resolution (m)
        if(.not.allocated(dx))            allocate(dx(nlat*nlon))
        !> \brief Grid cell resolution 2D (m)
        if(.not.allocated(dx_2d))         allocate(dx_2d(nlon,nlat))
        !> \brief Wind adjustment factor for wildfire applications
        if(.not.allocated(waf))           allocate(waf(nlat*nlon))
        !> \brief Wind adjustment factor for wildfire applications 2D
        if(.not.allocated(waf_2d))        allocate(waf_2d(nlon,nlat))
        !> \brief Flame height (m)
        if(.not.allocated(flameh))        allocate(flameh(nlat*nlon))
        !> \brief Flame height 2D (m)
        if(.not.allocated(flameh_2d))     allocate(flameh_2d(nlon,nlat))
    end if
    !> \}

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Diffusivity Profile Outputs
!-------------------------------------------------------------------------------

    !> \defgroup CanopyDiffusivity Canopy Diffusivity Arrays
    !! \brief Arrays for eddy diffusivity profiles within canopy
    !! \ingroup MemoryManagement
    !! \{
    if (ifcaneddy) then
        write(*,*)  'Canopy eddy Kz option selected'
        write(*,*)  '-------------------------------'
        !> \brief Eddy diffusivity profile within canopy (m² s⁻¹)
        if(.not.allocated(Kz))            allocate(Kz(nlat*nlon,modlays))
        !> \brief Eddy diffusivity profile within canopy 3D (m² s⁻¹)
        if(.not.allocated(Kz_3d))         allocate(Kz_3d(nlon,nlat,modlays))
    end if
    !> \}


!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Photolysis Correction Factor Outputs
!-------------------------------------------------------------------------------

    !> \defgroup CanopyPhotolysis Canopy Photolysis Arrays
    !! \brief Arrays for photolysis rate correction factors within canopy
    !! \ingroup MemoryManagement
    !! \{
    if (ifcanphot) then
        write(*,*)  'Canopy photolysis option selected'
        write(*,*)  '-------------------------------'
        !> \brief Photolysis rate correction factor profile within canopy (dimensionless)
        if(.not.allocated(rjcf))            allocate(rjcf(nlat*nlon,modlays))
        !> \brief Photolysis rate correction factor profile within canopy 3D (dimensionless)
        if(.not.allocated(rjcf_3d))         allocate(rjcf_3d(nlon,nlat,modlays))
    end if
    !> \}

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Biogenic Emissions Outputs
!-------------------------------------------------------------------------------

    !> \defgroup CanopyBiogenicEmissions Canopy Biogenic Emissions Arrays
    !! \brief Arrays for biogenic emission profiles within canopy
    !! \ingroup MemoryManagement
    !! \{
    if (ifcanbio) then
        write(*,*)  'Canopy biogenic emissions option selected'
        write(*,*)  '-------------------------------'
        if (biospec_opt == 0 .or. biospec_opt == 1) then
            !> \brief Isoprene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_isop))         allocate(emi_isop(nlat*nlon,modlays))
            !> \brief Isoprene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_isop_3d))      allocate(emi_isop_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 2) then
            !> \brief Myrcene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_myrc))         allocate(emi_myrc(nlat*nlon,modlays))
            !> \brief Myrcene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_myrc_3d))      allocate(emi_myrc_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 3) then
            !> \brief Sabinene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_sabi))         allocate(emi_sabi(nlat*nlon,modlays))
            !> \brief Sabinene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_sabi_3d))      allocate(emi_sabi_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 4) then
            !> \brief Limonene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_limo))         allocate(emi_limo(nlat*nlon,modlays))
            !> \brief Limonene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_limo_3d))      allocate(emi_limo_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 5) then
            !> \brief 3-Carene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_care))         allocate(emi_care(nlat*nlon,modlays))
            !> \brief 3-Carene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_care_3d))      allocate(emi_care_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 6) then
            !> \brief Ocimene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_ocim))         allocate(emi_ocim(nlat*nlon,modlays))
            !> \brief Ocimene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_ocim_3d))      allocate(emi_ocim_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 7) then
            !> \brief Beta-pinene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_bpin))         allocate(emi_bpin(nlat*nlon,modlays))
            !> \brief Beta-pinene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_bpin_3d))      allocate(emi_bpin_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 8) then
            !> \brief Alpha-pinene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_apin))         allocate(emi_apin(nlat*nlon,modlays))
            !> \brief Alpha-pinene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_apin_3d))      allocate(emi_apin_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 9) then
            !> \brief Other monoterpenes emission rate (kg m-3 s-1)
            if(.not.allocated(emi_mono))         allocate(emi_mono(nlat*nlon,modlays))
            !> \brief Other monoterpenes emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_mono_3d))      allocate(emi_mono_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 10) then
            !> \brief Farnesene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_farn))         allocate(emi_farn(nlat*nlon,modlays))
            !> \brief Farnesene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_farn_3d))      allocate(emi_farn_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 11) then
            !> \brief Caryophyllene emission rate (kg m-3 s-1)
            if(.not.allocated(emi_cary))         allocate(emi_cary(nlat*nlon,modlays))
            !> \brief Caryophyllene emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_cary_3d))      allocate(emi_cary_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 12) then
            !> \brief Other sesquiterpenes emission rate (kg m-3 s-1)
            if(.not.allocated(emi_sesq))         allocate(emi_sesq(nlat*nlon,modlays))
            !> \brief Other sesquiterpenes emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_sesq_3d))      allocate(emi_sesq_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 13) then
            !> \brief 2-methyl-3-buten-2-ol emission rate (kg m-3 s-1)
            if(.not.allocated(emi_mbol))         allocate(emi_mbol(nlat*nlon,modlays))
            !> \brief 2-methyl-3-buten-2-ol emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_mbol_3d))      allocate(emi_mbol_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 14) then
            !> \brief Methanol emission rate (kg m-3 s-1)
            if(.not.allocated(emi_meth))         allocate(emi_meth(nlat*nlon,modlays))
            !> \brief Methanol emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_meth_3d))      allocate(emi_meth_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 15) then
            !> \brief Acetone emission rate (kg m-3 s-1)
            if(.not.allocated(emi_acet))         allocate(emi_acet(nlat*nlon,modlays))
            !> \brief Acetone emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_acet_3d))      allocate(emi_acet_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 16) then
            !> \brief Carbon monoxide emission rate (kg m-3 s-1)
            if(.not.allocated(emi_co))           allocate(emi_co(nlat*nlon,modlays))
            !> \brief Carbon monoxide emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_co_3d))        allocate(emi_co_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 17) then
            !> \brief Bidirectional VOC emission rate (kg m-3 s-1)
            if(.not.allocated(emi_bvoc))         allocate(emi_bvoc(nlat*nlon,modlays))
            !> \brief Bidirectional VOC emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_bvoc_3d))      allocate(emi_bvoc_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 18) then
            !> \brief Stress VOC emission rate (kg m-3 s-1)
            if(.not.allocated(emi_svoc))         allocate(emi_svoc(nlat*nlon,modlays))
            !> \brief Stress VOC emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_svoc_3d))      allocate(emi_svoc_3d(nlon,nlat,modlays))
        end if
        if (biospec_opt == 0 .or. biospec_opt == 19) then
            !> \brief Other VOC emission rate (kg m-3 s-1)
            if(.not.allocated(emi_ovoc))         allocate(emi_ovoc(nlat*nlon,modlays))
            !> \brief Other VOC emission rate 3D (kg m-3 s-1)
            if(.not.allocated(emi_ovoc_3d))      allocate(emi_ovoc_3d(nlon,nlat,modlays))
        end if
    end if
    !> \}

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Gas Dry Deposition Outputs
!-------------------------------------------------------------------------------

    !> \defgroup CanopyDryDeposition Canopy Gas Dry Deposition Arrays
    !! \brief Arrays for gas dry deposition velocities within canopy
    !! \ingroup MemoryManagement
    !! \{
    if (ifcanddepgas) then
        write(*,*)  'Canopy gas dry deposition option selected'
        write(*,*)  '-------------------------------'
        if (chemmechgas_opt == 0) then !RACM2 --> 31 species
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 1) then
                !> \brief Nitric oxide dry deposition velocity (cm s⁻¹)
                if(.not.allocated(ddep_no))              allocate(ddep_no(nlat*nlon,modlays))
                !> \brief Nitric oxide dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_no_3d))           allocate(ddep_no_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 2) then
                !> \brief Nitrogen dioxide dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_no2))              allocate(ddep_no2(nlat*nlon,modlays))
                !> \brief Nitrogen dioxide dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_no2_3d))           allocate(ddep_no2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 3) then
                !> \brief Ozone dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_o3))              allocate(ddep_o3(nlat*nlon,modlays))
                !> \brief Ozone dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_o3_3d))           allocate(ddep_o3_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 4) then
                !> \brief Nitrous acid dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_hono))              allocate(ddep_hono(nlat*nlon,modlays))
                !> \brief Nitrous acid dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_hono_3d))           allocate(ddep_hono_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 5) then
                !> \brief Peroxynitric acid dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_hno4))              allocate(ddep_hno4(nlat*nlon,modlays))
                !> \brief Peroxynitric acid dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_hno4_3d))           allocate(ddep_hno4_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 6) then
                !> \brief Nitric acid dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_hno3))              allocate(ddep_hno3(nlat*nlon,modlays))
                !> \brief Nitric acid dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_hno3_3d))           allocate(ddep_hno3_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 7) then
                !> \brief Dinitrogen pentoxide dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_n2o5))              allocate(ddep_n2o5(nlat*nlon,modlays))
                !> \brief Dinitrogen pentoxide dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_n2o5_3d))           allocate(ddep_n2o5_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 8) then
                !> \brief Carbon monoxide dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_co))              allocate(ddep_co(nlat*nlon,modlays))
                !> \brief Carbon monoxide dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_co_3d))           allocate(ddep_co_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 9) then
                !> \brief Hydrogen peroxide dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_h2o2))              allocate(ddep_h2o2(nlat*nlon,modlays))
                !> \brief Hydrogen peroxide dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_h2o2_3d))           allocate(ddep_h2o2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 10) then
                !> \brief Methane dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_ch4))              allocate(ddep_ch4(nlat*nlon,modlays))
                !> \brief Methane dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_ch4_3d))           allocate(ddep_ch4_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 11) then
                !> \brief Methyl peroxy radical dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_mo2))              allocate(ddep_mo2(nlat*nlon,modlays))
                !> \brief Methyl peroxy radical dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_mo2_3d))           allocate(ddep_mo2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 12) then
                !> \brief Methyl hydroperoxide dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_op1))              allocate(ddep_op1(nlat*nlon,modlays))
                !> \brief Methyl hydroperoxide dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_op1_3d))           allocate(ddep_op1_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 13) then
                !> \brief Methanol dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_moh))              allocate(ddep_moh(nlat*nlon,modlays))
                !> \brief Methanol dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_moh_3d))           allocate(ddep_moh_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 14) then
                !> \brief Nitrate radical dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_no3))              allocate(ddep_no3(nlat*nlon,modlays))
                !> \brief Nitrate radical dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_no3_3d))           allocate(ddep_no3_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 15) then
                !> \brief Oxygen atom triplet dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_o3p))              allocate(ddep_o3p(nlat*nlon,modlays))
                !> \brief Oxygen atom triplet dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_o3p_3d))           allocate(ddep_o3p_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 16) then
                !> \brief Oxygen atom singlet dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_o1d))              allocate(ddep_o1d(nlat*nlon,modlays))
                !> \brief Oxygen atom singlet dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_o1d_3d))           allocate(ddep_o1d_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 17) then
                !> \brief Hydroxyl radical dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_ho))              allocate(ddep_ho(nlat*nlon,modlays))
                !> \brief Hydroxyl radical dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_ho_3d))           allocate(ddep_ho_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 18) then
                !> \brief Hydroperoxy radical dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_ho2))              allocate(ddep_ho2(nlat*nlon,modlays))
                !> \brief Hydroperoxy radical dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_ho2_3d))           allocate(ddep_ho2_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 19) then
                !> \brief Formic acid dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_ora1))              allocate(ddep_ora1(nlat*nlon,modlays))
                !> \brief Formic acid dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_ora1_3d))           allocate(ddep_ora1_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 20) then
                !> \brief Hydroxyacetone dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_hac))              allocate(ddep_hac(nlat*nlon,modlays))
                !> \brief Hydroxyacetone dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_hac_3d))           allocate(ddep_hac_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 21) then
                !> \brief Peroxyacetic acid dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_paa))              allocate(ddep_paa(nlat*nlon,modlays))
                !> \brief Peroxyacetic acid dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_paa_3d))           allocate(ddep_paa_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 22) then
                !> \brief Dihydroxymethyl butenone dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_dhmob))              allocate(ddep_dhmob(nlat*nlon,modlays))
                !> \brief Dihydroxymethyl butenone dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_dhmob_3d))           allocate(ddep_dhmob_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 23) then
                !> \brief Hydroxyperoxy aldehyde dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_hpald))              allocate(ddep_hpald(nlat*nlon,modlays))
                !> \brief Hydroxyperoxy aldehyde dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_hpald_3d))           allocate(ddep_hpald_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 24) then
                !> \brief Isoprene hydroxy hydroperoxide dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_ishp))              allocate(ddep_ishp(nlat*nlon,modlays))
                !> \brief Isoprene hydroxy hydroperoxide dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_ishp_3d))           allocate(ddep_ishp_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 25) then
                !> \brief Isoprene epoxydiol dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_iepox))              allocate(ddep_iepox(nlat*nlon,modlays))
                !> \brief Isoprene epoxydiol dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_iepox_3d))           allocate(ddep_iepox_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 26) then
                !> \brief Propanone nitrate dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_propnn))              allocate(ddep_propnn(nlat*nlon,modlays))
                !> \brief Propanone nitrate dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_propnn_3d))           allocate(ddep_propnn_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 27) then
                !> \brief Isoprene nitrate beta dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_isopnb))              allocate(ddep_isopnb(nlat*nlon,modlays))
                !> \brief Isoprene nitrate beta dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_isopnb_3d))           allocate(ddep_isopnb_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 28) then
                !> \brief Isoprene nitrate delta dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_isopnd))              allocate(ddep_isopnd(nlat*nlon,modlays))
                !> \brief Isoprene nitrate delta dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_isopnd_3d))           allocate(ddep_isopnd_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 29) then
                !> \brief Methacrolein nitrate dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_macrn))              allocate(ddep_macrn(nlat*nlon,modlays))
                !> \brief Methacrolein nitrate dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_macrn_3d))           allocate(ddep_macrn_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 30) then
                !> \brief Methyl vinyl ketone nitrate dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_mvkn))              allocate(ddep_mvkn(nlat*nlon,modlays))
                !> \brief Methyl vinyl ketone nitrate dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_mvkn_3d))           allocate(ddep_mvkn_3d(nlon,nlat,modlays))
            end if
            if (ddepspecgas_opt == 0 .or. ddepspecgas_opt == 31) then
                !> \brief Isoprene nitrate peroxy dry deposition velocity (cm s-1)
                if(.not.allocated(ddep_isnp))              allocate(ddep_isnp(nlat*nlon,modlays))
                !> \brief Isoprene nitrate peroxy dry deposition velocity 3D (cm s-1)
                if(.not.allocated(ddep_isnp_3d))           allocate(ddep_isnp_3d(nlon,nlat,modlays))
            end if
        else
            write(*,*)  'Wrong chemical mechanism option of ', chemmechgas_opt, ' in namelist...exiting'
            write(*,*)  'Set chemmechgas_opt to only 0 (RACM2) for now'
            call exit(2)
        end if
    end if
    !> \}

END SUBROUTINE canopy_alloc

!> \}
