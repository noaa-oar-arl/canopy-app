!> \file canopy_soilno_mod.F90
!! \brief BDSNP Soil NO Emissions Module
!! \details Module for computing soil NO emissions using the Berkeley-Dalhousie
!!          Soil NO Parameterization (BDSNP) model with canopy reduction factor.
!!
!!          The fertilizer N contribution currently uses a steady-state
!!          approximation: fert_flux = Ninput * fert_frac * unit_conversion.
!!          This is the steady-state limit of the Hudman et al. (2012) Eq. 5
!!          dynamic N mass-balance reservoir for constant annual input.
!!
!! \todo    **Seasonal N distribution & dynamic reservoir (future task)**:
!!          The current NPKGRIDS input is a single annual total with no temporal
!!          variation, so the steady-state formula is exact.  To capture transient
!!          fertilizer dynamics (post-application decay, growing-season shutoff),
!!          the following enhancements are needed:
!!          1. Add seasonal distribution of the NPKGRIDS annual N total (e.g.,
!!             using MODIS EVI phenology or a crop-calendar growing-season window)
!!             so that the application rate F varies month-to-month.
!!          2. Implement the full Hudman et al. (2012) Eq. 5 mass-balance:
!!             N_avail(t) = N_avail(t-dt)*exp(-dt/tau) + F*tau*(1-exp(-dt/tau))
!!             with a persistent N_reservoir array tracked across timesteps.
!!          3. Add namelist parameters: fert_tau (decay lifetime, default 4 months)
!!             and fert_runoff (runoff fraction, default 0.4).
!!
!! \author P. C. Campbell (Initial version, Oct 2025)
!! \author Quazi Rasool (CIRES/NOAA CSL) (Syntax fix, wiring, output integration, Feb 2026)


module canopy_soilno_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: compute_soil_no_emissions

contains

    subroutine compute_soil_no_emissions(Tsoil, Wsoil, Ninput, LAI, CRF_in, NO_flux, &
        VTYPE, LU_OPT, PRATE, WILT, crf_opt, fert_frac, LAD, ZK, FCH, MODLAYS, &
        TEMPA, PRESSA, RELHUMA, UBAR, FSUN, PPFD_SUN, PPFD_SHADE, &
        SRAD, D_H, HREF, UBZREF, TMPSURF, TMP2M, HCM, &
        CHEMMECHGAS_OPT, CHEMMECHGAS_TOT, RAMIN)
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
        ! Additional optional args for crf_opt=2,3 (BDSNP deposition-based CRF):
        !   TEMPA      - in-canopy temperature profile (K)
        !   PRESSA     - in-canopy pressure profile (mb)
        !   RELHUMA    - in-canopy relative humidity profile (%)
        !   UBAR       - in-canopy wind speed profile (m/s)
        !   FSUN       - sunlit fraction profile
        !   PPFD_SUN   - sunlit PPFD profile (umol/m2/s)
        !   PPFD_SHADE - shaded PPFD profile (umol/m2/s)
        !   SRAD       - solar radiation (W/m2)
        !   D_H        - displacement height / canopy height
        !   HREF       - reference height above canopy (m)
        !   UBZREF     - wind speed at reference height (m/s)
        !   TMPSURF    - surface temperature (K)
        !   TMP2M      - 2m temperature (K)
        !   HCM        - canopy height (m)
        !   CHEMMECHGAS_OPT - chemical mechanism option
        !   CHEMMECHGAS_TOT - total gas species count
        !   RAMIN      - minimum aerodynamic resistance (s/m)
        use canopy_utils_mod, only: MolecDiff, rs_zhang_gas, EffHenrysLawCoeff, &
            ReactivityParam, rbl, rcl, rml, rav, CalcRiB

        real(real64), intent(in)  :: Tsoil, Wsoil, Ninput, LAI, CRF_in, PRATE, WILT
        real(real64), intent(in)  :: fert_frac
        integer, intent(in)       :: crf_opt, MODLAYS
        real(real64), intent(in), optional :: LAD(:), ZK(:), FCH
        ! Optional args for crf_opt=2,3 (BDSNP deposition-based CRF)
        real(real64), intent(in), optional :: TEMPA(:), PRESSA(:), RELHUMA(:)
        real(real64), intent(in), optional :: UBAR(:), FSUN(:), PPFD_SUN(:), PPFD_SHADE(:)
        real(real64), intent(in), optional :: SRAD, D_H, HREF, UBZREF, TMPSURF, TMP2M, HCM
        real(real64), intent(in), optional :: RAMIN
        integer, intent(in), optional :: CHEMMECHGAS_OPT, CHEMMECHGAS_TOT
        integer, intent(in)       :: VTYPE, LU_OPT
        real(real64), intent(out) :: NO_flux

        ! BDSNP/MEGANv3.2 parameters
        real(real64), parameter :: Q10 = 2.0d0
        real(real64), parameter :: T_REF = 303.15d0 ! 30C
        real(real64), parameter :: WFPS_OPT = 0.20d0
        real(real64), parameter :: POROSITY = 0.5d0
        ! Unit conversion: kg-N/ha/yr -> ng-N/m2/s
        ! 1 kg = 1e12 ng, 1 ha = 1e4 m2, 1 yr = 3.1536e7 s
        ! => 1 kg-N/ha/yr = 1e12 / (1e4 * 3.1536e7) = 3.1710 ng-N/m2/s
        real(real64), parameter :: KG_HA_YR_TO_NG_M2_S = 3.1710d0
        ! Species index for NO2 in RACM2
        integer, parameter :: NO2_INDEX = 2

        real(real64) :: base_flux, fert_flux, T_factor, WFPS, WFPS_factor, pulse_factor
        real(real64) :: biome_emission, crf_val
        ! Variables for vertically resolved CRF (case 1)
        integer :: k, nlay
        real(real64) :: dz, attenuation
        ! Variables for BDSNP deposition-based CRF (case 2, 3)
        real(real64) :: RiB_loc, Ra_loc, hstarl, f01, mdiffl_k
        real(real64) :: rs_k, rb_k, rc_k, rm_k, rnum, rden, rlx, vd_k
        real(real64) :: ppfd_k, cum_lai, vd_bar, sum_vd, n_vd
        real(real64) :: lad_dz
        logical :: have_dep_args

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
        ! Ninput (kg-N/ha) * fert_frac * unit conversion
        ! fert_frac: 0.01 = 1% (Steinkamp & Lawrence 2011)
        !            0.025 = 2.5% (Hudman et al. 2012)
        !
        ! NOTE: This is the steady-state solution of Hudman et al. (2012) Eq. 5:
        !   N_avail(t) = N_avail(0)*exp(-t/tau) + F*tau*(1 - exp(-t/tau))
        ! For constant annual input F, N_avail -> F*tau as t -> infinity.
        ! Emission = (N_avail / tau) * fert_frac = F * fert_frac, which is
        ! exactly what this line computes.  The full reservoir with exponential
        ! decay (tau ~ 4 months) only adds value when F varies in time (e.g.,
        ! seasonal fertilizer application with growing-season shutoff).
        ! Current NPKGRIDS input is a single annual total with no temporal
        ! variation, so the reservoir would remain at steady state and produce
        ! identical results.
        ! TODO: When seasonal N distribution is added, implement full Eq. 5
        !       mass-balance with persistent N_reservoir array and decay lifetime.
        fert_flux = Ninput * fert_frac * KG_HA_YR_TO_NG_M2_S

        ! Base BDSNP soil NO emission (ng N m-2 s-1)
        ! = biome background (modulated by T, moisture, pulse) + fertilizer N
        base_flux = biome_emission * T_factor * WFPS_factor * pulse_factor + fert_flux


        select case (crf_opt)
          case (0)
            ! Default bulk CRF: Beer's law with k=0.5
            crf_val = exp(-0.5d0 * LAI)

          case (1)
            ! LAD-resolved Beer's law CRF using vertical LAD profile
            if (present(LAD) .and. present(ZK) .and. present(FCH)) then
                attenuation = 0.0d0
                nlay = min(MODLAYS, size(LAD))
                dz = ZK(2) - ZK(1)
                n_vd = 0.0d0
                do k = 1, nlay
                    if (ZK(k) > 0.0d0 .and. ZK(k) <= FCH) then
                        attenuation = attenuation + exp(-LAD(k)*dz)
                        n_vd = n_vd + 1.0d0
                    end if
                end do
                if (n_vd > 0.0d0) then
                    crf_val = attenuation / n_vd
                else
                    crf_val = exp(-0.5d0 * LAI)
                end if
                crf_val = max(0.1d0, min(1.0d0, crf_val))
            else
                crf_val = exp(-0.5d0 * LAI)
            end if

          case (2)
            ! BDSNP deposition-based CRF using bulk LAI
            ! CRF = exp(-Vd_NO2 * LAI / Vent)
            ! Physically: fraction of soil NO escaping canopy after
            ! oxidation to NO2 and subsequent leaf uptake.
            have_dep_args = present(TEMPA) .and. present(PRESSA) .and. &
                present(RELHUMA) .and. present(UBAR) .and. present(FSUN) .and. &
                present(PPFD_SUN) .and. present(PPFD_SHADE) .and. &
                present(SRAD) .and. present(D_H) .and. present(HREF) .and. &
                present(UBZREF) .and. present(TMPSURF) .and. present(TMP2M) .and. &
                present(HCM) .and. present(ZK) .and. present(FCH) .and. &
                present(CHEMMECHGAS_OPT) .and. present(CHEMMECHGAS_TOT) .and. &
                present(RAMIN)

            if (have_dep_args .and. HCM > 0.0d0) then
                RiB_loc = CalcRiB(TMP2M, TMPSURF, UBZREF*100.0d0, &
                    D_H*HCM*100.0d0, (HCM+HREF)*100.0d0)
                Ra_loc = rav(UBZREF*100.0d0, (HCM+HREF)*100.0d0, &
                    D_H*HCM*100.0d0, HCM*100.0d0, RiB_loc)
                Ra_loc = max(RAMIN/100.0d0, Ra_loc)

                hstarl = EffHenrysLawCoeff(CHEMMECHGAS_OPT, CHEMMECHGAS_TOT, NO2_INDEX)
                f01    = ReactivityParam(CHEMMECHGAS_OPT, CHEMMECHGAS_TOT, NO2_INDEX)

                sum_vd = 0.0d0
                n_vd   = 0.0d0
                nlay = min(MODLAYS, size(ZK))
                do k = 1, nlay
                    if (ZK(k) > 0.0d0 .and. ZK(k) <= FCH) then
                        mdiffl_k = MolecDiff(CHEMMECHGAS_OPT, CHEMMECHGAS_TOT, &
                            NO2_INDEX, real(TEMPA(k),8), real(PRESSA(k),8))
                        ppfd_k = PPFD_SUN(k)*FSUN(k) + PPFD_SHADE(k)*(1.0d0-FSUN(k))
                        rs_k = rs_zhang_gas(mdiffl_k, real(TEMPA(k),8), &
                            real(PRESSA(k),8), ppfd_k, SRAD, real(RELHUMA(k),8))
                        rb_k = rbl(mdiffl_k, UBAR(k)*100.0d0)
                        rc_k = rcl(hstarl, f01)
                        rm_k = rml(hstarl, f01)
                        rnum = rc_k * (rs_k + rm_k)
                        rden = rc_k + 2.0d0 * (rs_k + rm_k)
                        rlx  = rb_k + (rnum/rden) + Ra_loc
                        vd_k = 1.0d0 / rlx
                        sum_vd = sum_vd + vd_k
                        n_vd   = n_vd + 1.0d0
                    end if
                end do

                if (n_vd > 0.0d0) then
                    vd_bar = sum_vd / n_vd
                    crf_val = exp(-vd_bar * LAI / max(1.0d0, UBZREF*10.0d0))
                else
                    crf_val = exp(-0.5d0 * LAI)
                end if
                crf_val = max(0.1d0, min(1.0d0, crf_val))
            else
                crf_val = exp(-0.5d0 * LAI)
            end if

          case (3)
            ! BDSNP deposition-based CRF using LAD profile (layer-resolved)
            ! CRF = exp(-sum(vd_k * LAD_k * dz) / Vent)
            ! Like case 2 but uses per-layer LAD instead of bulk LAI
            have_dep_args = present(TEMPA) .and. present(PRESSA) .and. &
                present(RELHUMA) .and. present(UBAR) .and. present(FSUN) .and. &
                present(PPFD_SUN) .and. present(PPFD_SHADE) .and. &
                present(SRAD) .and. present(D_H) .and. present(HREF) .and. &
                present(UBZREF) .and. present(TMPSURF) .and. present(TMP2M) .and. &
                present(HCM) .and. present(ZK) .and. present(FCH) .and. &
                present(LAD) .and. &
                present(CHEMMECHGAS_OPT) .and. present(CHEMMECHGAS_TOT) .and. &
                present(RAMIN)

            if (have_dep_args .and. HCM > 0.0d0) then
                RiB_loc = CalcRiB(TMP2M, TMPSURF, UBZREF*100.0d0, &
                    D_H*HCM*100.0d0, (HCM+HREF)*100.0d0)
                Ra_loc = rav(UBZREF*100.0d0, (HCM+HREF)*100.0d0, &
                    D_H*HCM*100.0d0, HCM*100.0d0, RiB_loc)
                Ra_loc = max(RAMIN/100.0d0, Ra_loc)

                hstarl = EffHenrysLawCoeff(CHEMMECHGAS_OPT, CHEMMECHGAS_TOT, NO2_INDEX)
                f01    = ReactivityParam(CHEMMECHGAS_OPT, CHEMMECHGAS_TOT, NO2_INDEX)

                cum_lai = 0.0d0
                nlay = min(MODLAYS, size(ZK))
                dz = ZK(2) - ZK(1)
                do k = 1, nlay
                    if (ZK(k) > 0.0d0 .and. ZK(k) <= FCH) then
                        mdiffl_k = MolecDiff(CHEMMECHGAS_OPT, CHEMMECHGAS_TOT, &
                            NO2_INDEX, real(TEMPA(k),8), real(PRESSA(k),8))
                        ppfd_k = PPFD_SUN(k)*FSUN(k) + PPFD_SHADE(k)*(1.0d0-FSUN(k))
                        rs_k = rs_zhang_gas(mdiffl_k, real(TEMPA(k),8), &
                            real(PRESSA(k),8), ppfd_k, SRAD, real(RELHUMA(k),8))
                        rb_k = rbl(mdiffl_k, UBAR(k)*100.0d0)
                        rc_k = rcl(hstarl, f01)
                        rm_k = rml(hstarl, f01)
                        rnum = rc_k * (rs_k + rm_k)
                        rden = rc_k + 2.0d0 * (rs_k + rm_k)
                        rlx  = rb_k + (rnum/rden) + Ra_loc
                        vd_k = 1.0d0 / rlx
                        lad_dz = LAD(k) * dz
                        cum_lai = cum_lai + vd_k * lad_dz
                    end if
                end do

                if (cum_lai > 0.0d0) then
                    crf_val = exp(-cum_lai / max(1.0d0, UBZREF*10.0d0))
                else
                    crf_val = exp(-0.5d0 * LAI)
                end if
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
