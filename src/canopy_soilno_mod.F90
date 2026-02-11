!> \file canopy_soilno_mod.F90
!! \brief BDSNP Soil NO Emissions Module
!! \details Module for computing soil NO emissions using the Berkeley-Dalhousie
!!          Soil NO Parameterization (BDSNP) model with canopy reduction factor.
!! \author P. C. Campbell (Initial version, Oct 2025)
!! \author Quazi Rasool (CIRES/NOAA CSL) (Syntax fix, wiring, output integration, Feb 2026)


module canopy_soilno_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: compute_soil_no_emissions

contains

    subroutine compute_soil_no_emissions(Tsoil, Wsoil, Ninput, LAI, CRF_in, NO_flux, &
        VTYPE, LU_OPT, PRATE, WILT, crf_opt, fert_frac, LAD, ZK, FCH, MODLAYS)
        ! Compute soil NO emissions using BDSNP model and robust CRF
        ! Arguments:
        !   Tsoil      - soil temperature (K)
        !   Wsoil      - soil water content (m3/m3)
        !   Ninput     - nitrogen input (kg N/ha)
        !   LAI        - leaf area index (m2/m2)
        !   CRF_in     - user override for canopy reduction factor (0=auto)
        !   NO_flux    - output: soil NO flux (ng N m-2 s-1)
        !   VTYPE      - vegetation type (integer)
        !   LU_OPT     - land use option (integer)
        !   PRATE      - precipitation rate (mm/hr)
        !   WILT       - wilting point (m3/m3)
        !   fert_frac  - fraction of applied N emitted as NO
        !                0.01 = 1% (Steinkamp & Lawrence 2011)
        !                0.025 = 2.5% (Hudman et al. 2012)
        real(real64), intent(in)  :: Tsoil, Wsoil, Ninput, LAI, CRF_in, PRATE, WILT
        real(real64), intent(in)  :: fert_frac
        integer, intent(in)       :: crf_opt, MODLAYS
        real(real64), intent(in), optional :: LAD(:), ZK(:), FCH
        integer, intent(in)       :: VTYPE, LU_OPT
        real(real64), intent(out) :: NO_flux

        ! BDSNP/MEGANv3.2 parameters
        real(real64), parameter :: Q10 = 2.0d0
        real(real64), parameter :: T_REF = 303.15d0 ! 30C
        real(real64), parameter :: WFPS_OPT = 0.20d0
        real(real64), parameter :: POROSITY = 0.5d0
        ! Unit conversion: kg-N/ha/yr -> ng-N/m2/s
        ! 1 kg = 1e12 ng, 1 ha = 1e4 m2, 1 yr = 3.1536e7 s
        ! => 1 kg-N/ha/yr = 1e12 / (1e4 * 3.1536e7) = 3.1709d-3 ng-N/m2/s
        real(real64), parameter :: KG_HA_YR_TO_NG_M2_S = 3.1709d-3

        real(real64) :: base_flux, fert_flux, T_factor, WFPS, WFPS_factor, pulse_factor
        real(real64) :: biome_emission, crf_val
        ! Variables for vertically resolved CRF (case 2)
        integer :: k, nlay
        real(real64) :: dz, sum_lad, attenuation

        ! Calculate water-filled pore space (WFPS)
        WFPS = min(1.0d0, Wsoil / POROSITY)

        ! Get biome-specific emission factor (ng N m-2 s-1)
        biome_emission = get_biome_emission(VTYPE, LU_OPT)

        ! Temperature response
        T_factor = Q10 ** ((Tsoil - T_REF) / 10.0d0)

        ! Soil moisture response (Guenther et al. 2006, MEGANv3.2)
        if (WFPS <= WFPS_OPT) then
            WFPS_factor = WFPS / WFPS_OPT
        else
            WFPS_factor = 1.0d0 - ((WFPS - WFPS_OPT) / (1.0d0 - WFPS_OPT))
        end if
        WFPS_factor = max(0.0d0, min(1.0d0, WFPS_factor))

        ! Pulsing factor (rain pulse, simplified)
        pulse_factor = 1.0d0
        if (PRATE > 0.0d0 .and. WFPS < 0.3d0) then
            pulse_factor = 10.7d0 * WFPS + 1.0d0
        end if

        ! Fertilizer N contribution (ng N m-2 s-1)
        ! Ninput (kg-N/ha/yr) * fert_frac * unit conversion
        ! fert_frac: 0.01 = 1% (Steinkamp & Lawrence 2011)
        !            0.025 = 2.5% (Hudman et al. 2012)
        fert_flux = Ninput * fert_frac * KG_HA_YR_TO_NG_M2_S

        ! Base BDSNP soil NO emission (ng N m-2 s-1)
        ! = biome background (modulated by T, moisture, pulse) + fertilizer N
        base_flux = biome_emission * T_factor * WFPS_factor * pulse_factor + fert_flux


        select case (crf_opt)
          case (0)
            ! Default CRF
            crf_val = exp(-0.5d0 * LAI)
          case (1)
            ! Robust CRF (same as default for now, can be refined)
            crf_val = exp(-0.5d0 * LAI)
          case (2)
            ! Vertically resolved CRF using LAD profile (if provided)
            if (present(LAD) .and. present(ZK) .and. present(FCH)) then
                ! Integrate LAD over canopy height to get total attenuation
                sum_lad = 0.0d0
                attenuation = 0.0d0
                nlay = min(MODLAYS, size(LAD))
                do k = 1, nlay
                    if (ZK(k) <= FCH) then
                        dz = ZK(2)-ZK(1)
                        sum_lad = sum_lad + LAD(k)*dz
                        attenuation = attenuation + exp(-LAD(k)*dz)
                    end if
                end do
                ! CRF as mean attenuation across canopy
                crf_val = attenuation / nlay
                crf_val = max(0.1d0, min(1.0d0, crf_val))
            else
                crf_val = exp(-0.5d0 * LAI)
            end if
          case default
            crf_val = exp(-0.5d0 * LAI)
        end select
        if (CRF_in > 0.0d0) crf_val = CRF_in
        crf_val = max(0.1d0, min(1.0d0, crf_val))

        ! Final soil NO flux (ng N m-2 s-1)
        NO_flux = base_flux * crf_val
    end subroutine compute_soil_no_emissions

    function get_biome_emission(VTYPE, LU_OPT) result(emission_rate)
        ! Returns base soil NO emission rate for a specific biome/vegetation type (ng N m-2 s-1)
        integer, intent(in) :: VTYPE, LU_OPT
        real(real64) :: emission_rate
        ! Default values based on Hudman et al. (2012), HEMCO, MEGANv3.2
        emission_rate = 1.0d0
        if (LU_OPT == 0) then
            select case (VTYPE)
              case (1); emission_rate = 2.0d0   ! Evergreen Needleleaf Forest
              case (2); emission_rate = 2.5d0   ! Evergreen Broadleaf Forest
              case (3); emission_rate = 2.0d0   ! Deciduous Needleleaf Forest
              case (4); emission_rate = 2.5d0   ! Deciduous Broadleaf Forest
              case (5); emission_rate = 2.3d0   ! Mixed Forests
              case (6); emission_rate = 1.5d0   ! Closed Shrublands
              case (7); emission_rate = 1.0d0   ! Open Shrublands
              case (8); emission_rate = 3.0d0   ! Woody Savannas
              case (9); emission_rate = 4.0d0   ! Savannas
              case (10); emission_rate = 5.0d0  ! Grasslands
              case (11); emission_rate = 0.2d0  ! Permanent Wetlands
              case (12); emission_rate = 7.0d0  ! Croplands
              case (13); emission_rate = 0.5d0  ! Urban
              case (14); emission_rate = 6.0d0  ! Cropland/Natural Mosaic
              case (15); emission_rate = 0.0d0  ! Snow/Ice
              case (16); emission_rate = 0.2d0  ! Barren
              case default; emission_rate = 1.0d0
            end select
        else
            select case (VTYPE)
              case (1:5); emission_rate = 2.5d0
              case (6:7); emission_rate = 1.5d0
              case (8:9); emission_rate = 3.5d0
              case (10); emission_rate = 5.0d0
              case (12); emission_rate = 7.0d0
              case default; emission_rate = 1.0d0
            end select
        end if
    end function get_biome_emission

end module canopy_soilno_mod
