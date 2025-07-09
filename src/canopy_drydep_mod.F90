!> \file canopy_drydep_mod.F90
!> \brief Gas dry deposition calculations for canopy model
!> \details This module implements parameterized canopy gas dry deposition
!!          based on Zhang et al. (2003) algorithms for various surface types
!!          including vegetation, soil, snow, urban surfaces, and water.
!> \author P.C. Campbell
!> \date February 2025
!> \version 1.0

!> \defgroup drydep_group Gas Dry Deposition Calculations
!> \brief Routines for calculating gas dry deposition to various surface types
!> \details This group contains subroutines for calculating gas dry deposition
!!          velocities to vegetation, soil, snow, urban surfaces, and water
!!          based on Zhang et al. (2003) and other established parameterizations.
!> \{

module canopy_drydep_mod

    implicit none

contains

!> \brief Gas dry deposition to vegetation using Zhang et al. (2003) parameterization
!> \details Computes parameterized canopy gas dry deposition based on Zhang et al. (2003)
!!          algorithms as implemented in ACCESS (Saylor 2013). Calculates deposition
!!          velocities for gas species to vegetation surfaces within the canopy.
!> \param CHEMMECHGAS_OPT Chemical mechanism option selector
!> \param CHEMMECHGAS_TOT Total number of gas species in chemical mechanism
!> \param ZK Model heights (m)
!> \param FCH Canopy height (m)
!> \param TEMPA Ambient temperature profile in canopy (K)
!> \param PRESSA Ambient pressure profile in canopy (mb)
!> \param RELHUMA Ambient relative humidity profile in canopy (%)
!> \param FSUN Sunlit/shaded fraction from photolysis correction factor
!> \param PPFD_SUN PPFD for sunlit leaves (umol phot/m2 s)
!> \param PPFD_SHADE PPFD for shaded leaves (umol phot/m2 s)
!> \param UBAR Mean above/in-canopy wind speed (m/s)
!> \param SRAD Incoming solar irradiation top of canopy (W/m^2)
!> \param RA Aerodynamic resistance (s/cm)
!> \param DEP_IND Gas deposition species index (depends on gas mechanism)
!> \param DEP_OUT Output canopy layer gas dry deposition rate (cm/s) [output]
!> \author P.C. Campbell
!> \date February 2025
!> \note Based on Zhang et al. (2003) and Saylor (2013) ACCESS model
    SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        ZK, FCH, TEMPA, PRESSA, &
        RELHUMA, FSUN, PPFD_SUN, PPFD_SHADE, UBAR, &
        SRAD, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !< Constants for canopy models
        use canopy_utils_mod,  ONLY: MolecDiff,rs_zhang_gas,EffHenrysLawCoeff,& !< Utility functions
            ReactivityParam,rbl,rcl,rml

        !> \name Input Parameters
        !> \{
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    !< Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    !< Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: ZK(:)              !< Model heights (m)
        REAL(RK),    INTENT( IN )       :: FCH                !< Canopy height (m)
        REAL(RK),    INTENT( IN )       :: FSUN(:)            !< Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( IN )       :: PPFD_SUN(:)        !< PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: PPFD_SHADE(:)      !< PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( IN )       :: TEMPA(:)           !< Ambient Temperature profile in canopy (K)
        REAL(RK),    INTENT( IN )       :: PRESSA(:)          !< Ambient Pressure profile in canopy (mb)
        REAL(RK),    INTENT( IN )       :: RELHUMA(:)         !< Ambient Relative Humidity profile in canopy (%)
        REAL(RK),    INTENT( IN )       :: UBAR(:)            !< Mean above/in-canopy wind speed (m/s)
        REAL(RK),    INTENT( IN )       :: SRAD               !< Incoming solar irradiation top of canopy (W/m^2)
        REAL(RK),    INTENT( IN )       :: RA                 !< Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            !< Gas deposition species index (depends on gas mech, set in constants)
        !> \}

        !> \name Output Parameters
        !> \{
        REAL(RK),    INTENT( OUT )      :: DEP_OUT(:)         !< Output canopy layer gas dry deposition rate for each DEP_IND (cm/s)
        !> \}

        !> \name Local Variables
        !> \{
        REAL(RK) ::  PPFD(SIZE(ZK))                           !< PPFD ave sun and shade (umol/m2 s)
        REAL(RK) ::  mdiffl(SIZE(ZK))                         !< Molecular diffusivity for species l based on DEP_IND (cm2/s)
        REAL(RK) ::  rs(SIZE(ZK))                             !< Stomatal resistance for species l based on DEP_IND (s/cm)
        REAL(RK) ::  rb(SIZE(ZK))                             !< Boundary layer resistance for species l based on DEP_IND (s/cm)
        REAL(RK) ::  rc(SIZE(ZK))                             !< Cuticular resistance for species l based on DEP_IND (s/cm)
        REAL(RK) ::  rm(SIZE(ZK))                             !< Mesophyll resistance for species l based on DEP_IND (s/cm)
        REAL(RK) ::  hstarl                                   !< Effective Henry's law coefficient based on DEP_IND (M/atm)
        REAL(RK) ::  f01                                      !< Reactivity parameter based on DEP_IND (0-1)
        REAL(RK) :: rnum,rden,rlx,vdlx                        !< Working variables for resistance calculations
        INTEGER i                                             !< Loop index
        !> \}

        !> Average PPFD = sum sun and shade weighted by sunlit fraction
        PPFD = (PPFD_SUN*FSUN) + (PPFD_SHADE*(1.0-FSUN))

        !> Calculate molecular diffusivity (cm^2/s) and resistances (s/cm) of species l using from DEP_IND
        hstarl  = EffHenrysLawCoeff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)
        f01     = ReactivityParam(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)
        do i=1, SIZE(ZK)
            if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then           !< Above ground level and at/below canopy top
                mdiffl(i)  = MolecDiff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND,TEMPA(i),PRESSA(i))
                rs(i)      = rs_zhang_gas(mdiffl(i),TEMPA(i),PRESSA(i),PPFD(i),SRAD,RELHUMA(i)) !< Stomatal resistance (s/cm)
                rb(i)      = rbl(mdiffl(i), UBAR(i)*100.0_rk) !< Leaf boundary layer resistance (s/cm)
                rc(i)      = rcl(hstarl, f01)                              !< Leaf cuticular resistance (s/cm)
                rm(i)      = rml(hstarl, f01)                              !< Leaf mesophyll resistance (s/cm)
                rnum = rc(i) * (rs(i) + rm(i))
                rden = rc(i) + 2.0_rk * (rs(i) + rm(i))
                rlx   = rb(i) + (rnum/rden) + RA                           !< Surface+boundary+aerodynamic resistances (s/cm)
                vdlx  = 1.0_rk/rlx
                dep_out(i) = vdlx                                          !< Calculate deposition velocity (cm/s)
            else
                rb(i) = 0.0_rk
                rc(i) = 0.0_rk
                rm(i) = 0.0_rk
                rs(i) = 0.0_rk
                dep_out(i) = 0.0_rk
            endif
        end do

    END SUBROUTINE CANOPY_GAS_DRYDEP_ZHANG

!> \brief Gas dry deposition to soil surfaces
!> \details Calculates gas dry deposition velocity to soil surfaces based on
!!          molecular diffusivity and soil properties including soil type,
!!          depth, and moisture content.
!> \param CHEMMECHGAS_OPT Chemical mechanism option selector
!> \param CHEMMECHGAS_TOT Total number of gas species in chemical mechanism
!> \param TEMPSOIL Soil temperature in topsoil (K)
!> \param PRESSA Ambient pressure just above surface (mb)
!> \param UBAR Mean wind speed just above surface (m/s)
!> \param SOCAT Input soil category dataset used
!> \param SOTYP Input soil type integer associated with soil category
!> \param DSOIL Depth of topsoil (cm)
!> \param STHETA Volumetric soil water content in topsoil (m^3/m^3)
!> \param RA Aerodynamic resistance (s/cm)
!> \param DEP_IND Gas deposition species index (depends on gas mechanism)
!> \param DEP_OUT Output soil layer gas dry deposition rate (cm/s) [output]
!> \author P.C. Campbell
!> \date February 2025
    SUBROUTINE CANOPY_GAS_DRYDEP_SOIL( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        TEMPSOIL, PRESSA, UBAR, SOCAT, SOTYP, DSOIL, STHETA, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !< Constants for canopy models
        use canopy_utils_mod,  ONLY: MolecDiff,SoilResist,SoilRbg !< Utility functions

        !> \name Input Parameters
        !> \{
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    !< Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    !< Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: TEMPSOIL           !< Soil temperature in topsoil (K)
        REAL(RK),    INTENT( IN )       :: PRESSA             !< Ambient Pressure just above surface (mb)
        REAL(RK),    INTENT( IN )       :: UBAR               !< Mean wind speed just above surface (m/s)
        INTEGER,     INTENT( IN )       :: SOCAT              !< Input soil category dataset used
        INTEGER,     INTENT( IN )       :: SOTYP              !< Input soil type integer associated with soilcat
        REAL(RK),    INTENT( IN )       :: DSOIL              !< Depth of topsoil (cm)
        REAL(RK),    INTENT( IN )       :: STHETA             !< Volumetric soil water content in topsoil (m^3/m^3)
        REAL(RK),    INTENT( IN )       :: RA                 !< Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            !< Gas deposition species index (depends on gas mech)
        !> \}

        !> \name Output Parameters
        !> \{
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            !< Output soil layer gas dry deposition rate for each DEP_IND (cm/s)
        !> \}

        !> \name Local Variables
        !> \{
        real(rk)                        :: mdiffl             !< Molecular diffusivity (cm^2/s)
        real(rk)                        :: rsoill             !< Resistance to diffusion thru soil pore space for chemical species (s/cm)
        real(rk)                        :: rbg                !< Ground boundary layer resistance (s/cm)
        !> \}


        !> Use soil temperature and pressure
        mdiffl = MolecDiff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND,TEMPSOIL,PRESSA)

        !> Depends on soil type, depth, and moisture
        rsoill = SoilResist(mdiffl,SOCAT,SOTYP,DSOIL,STHETA)

        !> Convert wind to cm/s - Ground boundary layer resistance (s/cm)
        rbg = SoilRbg(UBAR*100.0_rk)
        !> Rbg is invariant to species not layers
        !! Must use model layer above surface as physically correct no-slip boundary condition is applied for wind speed at z = 0
        !> Deposition velocity to ground surface under canopy or outside of canopy, e.g., barren land (cm/s)
        DEP_OUT = 1.0_rk/(rbg+rsoill+RA)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_SOIL

!> \brief Gas dry deposition to snow surfaces
!> \details Calculates gas dry deposition velocity to snow surfaces based on
!!          reactivity relative to HNO3 following CMAQv5.3.1 formulation
!!          and Helmig et al. parameterizations.
!> \param CHEMMECHGAS_OPT Chemical mechanism option selector
!> \param CHEMMECHGAS_TOT Total number of gas species in chemical mechanism
!> \param UBAR Mean wind speed just above surface (m/s)
!> \param RA Aerodynamic resistance (s/cm)
!> \param DEP_IND Gas deposition species index (depends on gas mechanism)
!> \param DEP_OUT Output snow layer gas dry deposition rate (cm/s) [output]
!> \author P.C. Campbell
!> \date February 2025
!> \note Based on CMAQv5.3.1 formulation scaled to reactivity relative to HNO3
    SUBROUTINE CANOPY_GAS_DRYDEP_SNOW( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, UBAR, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !< Constants for canopy models
        use canopy_utils_mod,  ONLY: ReactivityParamHNO3, SoilRbg !< Utility functions

        !> \name Input Parameters
        !> \{
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    !< Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    !< Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: UBAR               !< Mean wind speed just above surface (m/s)
        REAL(RK),    INTENT( IN )       :: RA                 !< Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            !< Gas deposition species index (depends on gas mech)
        !> \}

        !> \name Output Parameters
        !> \{
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            !< Output soil layer gas dry deposition rate for each DEP_IND (cm/s)
        !> \}

        !> \name Local Variables and Parameters
        !> \{
        real(rk), parameter             :: ar_0   = 8.0       !< Used to scale other species to HNO3 (dimensionless)
        real(rk)                        :: ar_l               !< Reactivity denominator relative to HNO3 for each species (dimensionless)
        real(rk), parameter             :: rsnow0 = 100.0     !< Resistance to deposition to snow (s/cm) based on Helmig et al.
        real(rk)                        :: rsnowl             !< Resistance to diffusion thru snow space for chemical species (s/cm)
        real(rk)                        :: rbg                !< Ground boundary layer resistance (s/cm)
        !> \}

        ar_l = ReactivityParamHNO3(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        !> Based on CMAQv5.3.1 formulation scaled to reactivity relative to HNO3
        rsnowl = rsnow0 * (ar_0/ar_l)

        !> Convert wind to cm/s - Ground boundary layer resistance (s/cm)
        rbg = SoilRbg(UBAR*100.0_rk)
        !> Rbg is invariant to species not layers
        !! Must use model layer above surface as physically correct no-slip boundary condition is applied for wind speed at z = 0

        !> Deposition velocity to ground surface under canopy or outside of canopy, e.g., barren land (cm/s) when snow cover is present
        DEP_OUT = 1.0_rk/(rbg+rsnowl+RA)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_SNOW

!> \brief Gas dry deposition to urban surfaces
!> \details Calculates gas dry deposition velocity to urban surfaces based on
!!          reaction probability with building materials and Maxwell-Boltzmann
!!          gas kinetics following Shen and Gao (2018) parameterization.
!> \param CHEMMECHGAS_OPT Chemical mechanism option selector
!> \param CHEMMECHGAS_TOT Total number of gas species in chemical mechanism
!> \param UBAR Mean wind speed just above surface (m/s)
!> \param TEMP Mean temperature just above surface (K)
!> \param GAMMA_BUILD Reaction probability with building type (dimensionless, default 5.0e-5)
!> \param RA Aerodynamic resistance (s/cm)
!> \param DEP_IND Gas deposition species index (depends on gas mechanism)
!> \param DEP_OUT Output urban surface gas dry deposition rate (cm/s) [output]
!> \author P.C. Campbell
!> \date February 2025
!> \note Based on Shen and Gao (2018) parameterization
    SUBROUTINE CANOPY_GAS_DRYDEP_URBAN( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, UBAR, TEMP, GAMMA_BUILD, &
        RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk, pi, rgasuniv                       !< Constants for canopy models
        use canopy_utils_mod,  ONLY: MolarMassGas, SoilRbg !< Utility functions

        !> \name Input Parameters
        !> \{
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT    !< Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT    !< Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: UBAR               !< Mean wind speed just above surface (m/s)
        REAL(RK),    INTENT( IN )       :: TEMP               !< Mean temperature just above surface (K)
        REAL(RK),    INTENT( IN )       :: GAMMA_BUILD        !< Reaction probability with building type (dimensionless)
        !! Default NL is average of range in gamma from as low as 10−8 for
        !! glass and metal to 10−4 for activated carbon and brick.
        !! =5.0D-5.  Reference (Shen and Gao, 2018;
        !! https://doi.org/10.1016/j.buildenv.2018.02.046)
        REAL(RK),    INTENT( IN )       :: RA                 !< Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            !< Gas deposition species index (depends on gas mech)
        !> \}

        !> \name Output Parameters
        !> \{
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            !< Output soil layer gas dry deposition rate for each DEP_IND (cm/s)
        !> \}

        !> \name Local Variables
        !> \{
        real(rk)                        :: mmg_l              !< Molar mass for each gas species (kg/mol)
        real(rk)                        :: cave_l             !< Maxwell-Boltzmann average speed of gas distribution (m/s)
        real(rk)                        :: rurbanl            !< Resistance to diffusion thru snow space for chemical species (s/m)
        real(rk)                        :: rbg                !< Ground boundary layer resistance (s/cm)
        !> \}

        !> Get molar mass for each gas species (kg/mol)
        mmg_l = MolarMassGas(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        !> Based on Maxwell-Boltzmann distribution, average speed of gas distribution (cm/s)
        cave_l = sqrt((8.0_rk*rgasuniv*TEMP)/(pi*mmg_l))*100.0_rk

        !> Based on Shen and Gao (2018), Eq. (2):  https://doi.org/10.1016/j.buildenv.2018.02.046)
        rurbanl = 4.0_rk/(GAMMA_BUILD*cave_l) !< Already converted in to units of s/cm from cave_l

        !> Convert wind to cm/s - Ground boundary layer resistance (s/cm)
        rbg = SoilRbg(UBAR*100.0_rk)
        !> Rbg is invariant to species not layers
        !! Must use model layer above surface as physically correct no-slip boundary condition is applied for wind speed at z = 0

        !> Deposition velocity to urban surfaces (cm/s)
        DEP_OUT = 1.0_rk/(rbg+rurbanl+RA)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_URBAN

!> \brief Gas dry deposition to water surfaces
!> \details Calculates gas dry deposition velocity to water surfaces based on
!!          Henry's law coefficient, wet bulb temperature calculations, and
!!          diffusivity considerations following established parameterizations.
!> \param CHEMMECHGAS_OPT Chemical mechanism option selector
!> \param CHEMMECHGAS_TOT Total number of gas species in chemical mechanism
!> \param TEMP2 Mean temperature just above surface (K)
!> \param QV2 Mean mixing ratio just above surface (kg/kg)
!> \param USTAR Friction velocity at surface (m/s)
!> \param RA Aerodynamic resistance (s/cm)
!> \param DEP_IND Gas deposition species index (depends on gas mechanism)
!> \param DEP_OUT Output water surface gas dry deposition rate (cm/s) [output]
!> \author P.C. Campbell
!> \date February 2025
!> \note Based on Slinn et al. (1978), Fairall et al. (2007), and Fritschen and Gay (1979)
    SUBROUTINE CANOPY_GAS_DRYDEP_WATER( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, TEMP2, QV2, &
        USTAR, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk, cpd, lv0, dlvdt, stdtemp      !< Constants for canopy models
        use canopy_utils_mod,  ONLY: EffHenrysLawCoeff, LeBasMVGas, WaterRbw !< Utility functions

        !> \name Input Parameters
        !> \{
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_OPT     !< Select chemical mechanism
        INTEGER,     INTENT( IN )       :: CHEMMECHGAS_TOT     !< Select chemical mechanism gas species list
        REAL(RK),    INTENT( IN )       :: TEMP2               !< Mean temperature just above surface (K)
        REAL(RK),    INTENT( IN )       :: QV2                 !< Mean mixing ratio just above surface (kg/kg)
        REAL(RK),    INTENT( IN )       :: USTAR               !< Friction velocity at surface (m/s)
        REAL(RK),    INTENT( IN )       :: RA                 !< Aerodynamic resistance (s/cm)
        INTEGER,     INTENT( IN )       :: DEP_IND            !< Gas deposition species index (depends on gas mech)
        !> \}

        !> \name Output Parameters
        !> \{
        REAL(RK),    INTENT( OUT )      :: DEP_OUT            !< Output soil layer gas dry deposition rate for each DEP_IND (cm/s)
        !> \}

        !> \name Local Variables
        !> \{
        real(rk)                        :: hstarl             !< Effective Henry's law coefficient based on DEP_IND (M/atm)
        real(rk)                        :: ctemp2             !< Mean temperature just above surface (C)
        real(rk)                        :: lv                 !< Latent heat of vaporization (J/kg)
        real(rk)                        :: cp_air             !< Specific heat of moist air (J/kg-K)
        real(rk)                        :: tw                 !< Wet bulb temperature (K)
        real(rk)                        :: lebas_l            !< Le Bas molar volumes are from the Schroeder additive method (cm3/mol)
        real(rk)                        :: dw                 !< Diffusivity of water
        real(rk)                        :: dw25               !< Diffusivity of water at 298.15 K
        real(rk)                        :: kvisw              !< Kinematic viscosity of water (cm^2/s)
        real(rk)                        :: scw_pr_23          !< (scw/pr)**2/3
        real(rk)                        :: rbw                !< Water boundary layer resistance (s/cm)
        real(rk)                        :: rwaterl            !< Water surface resistance (s/cm)
        !> \}

        !> \name Physical Constants
        !> \{
        real(rk), Parameter             :: pr         = 0.709 !< Prandtl Number (dimensionless)
        real(rk), Parameter             :: rt25inK    = 1.0_rk/(stdtemp + 25.0_rk) !< 298.15K = 25C
        real(rk), Parameter             :: twothirds  = 2.0_rk / 3.0_rk !< Two thirds constant
        real(rk), Parameter             :: d3         = 1.38564e-2 !< Scaling parameter used to estimate the friction velocity in surface waters from
        !! the atmospheric friction velocity to a value following Slinn et al. (1978)
        !! and Fairall et al. (2007)
        !> \}

        hstarl  = EffHenrysLawCoeff(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        !> Calculate the water surface film temperature: wet bulb temperature.
        !! Wet bulb temperature based on eqn in Fritschen and Gay (1979).
        ctemp2 = TEMP2 - stdtemp                                      !< Convert to C
        lv     = lv0 - dlvdt * ctemp2                                 !< Latent heat (J/kg)
        cp_air = cpd * ( 1.0_rk + 0.84_rk * QV2 )                     !< Specific heat (J/kg-K)
        tw     = ( ( 4.71e4 * cp_air / lv ) - 0.870_rk ) + stdtemp    !< Wet bulb temp (K)

        !> Make Henry's Law constant non-dimensional.
        hstarl  = hstarl * 0.08205_rk * tw

        !> Get Le Bas molar volumes for each gas species (cm3/mol)
        lebas_l = LeBasMVGas(CHEMMECHGAS_OPT,CHEMMECHGAS_TOT,DEP_IND)

        !> From Hayduk and Laudie parameterization
        dw25 = 13.26e-5 / ( 0.8904_rk**1.14_rk * lebas_l**0.589_rk )
        kvisw = 0.017_rk * EXP( -0.025_rk * ( tw - stdtemp ) )
        dw    = dw25 * ( tw * rt25inK ) * ( 0.009025_rk / kvisw )
        scw_pr_23 = ( ( kvisw / dw ) / pr ) ** twothirds

        !> Resistance to water surface (s/cm)
        rwaterl = scw_pr_23 / ( hstarl * d3 * USTAR*100.0_rk )

        !> Convert ustar to cm/s and scale by d3 - Water boundary layer resistance (s/cm)
        rbw = WaterRbw(d3*USTAR*100.0_rk)
        !> Rbw is invariant to species not layers

        !> Deposition velocity to water surfaces (cm/s)
        DEP_OUT = 1.0_rk/(rbw+rwaterl+RA)

        return
    END SUBROUTINE CANOPY_GAS_DRYDEP_WATER

!> \}

end module canopy_drydep_mod
