module canopy_bioemi_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_BIO( ZK, FCLAI, FCH, LAI, CLU, COSZEN, SFCRAD, &
        TEMP2, LU_OPT, VTYPE, EMI_IND, EMI_OUT)

!-----------------------------------------------------------------------

! Description:
!     computes parameterized canopy biogenic emissions

! Preconditions:
!     in-canopy FCLAI, model LAI, surface/skin temperature.

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Clifton et al. (2022) algorithms
! Citation:
! Clifton, O. E. et al. (2022). Large eddy simulation for investigating
! coupled forest canopy and turbulence influences on atmospheric chemistry.
! Journal of Advances in Modeling Earth Systems, 14, e2022MS003078.
! https://doi.org/10.1029/2022MS003078
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Jan 2023 P.C. Campbell: Initial canopy isoprene only version
!     Feb 2023 P.C. Campbell: Modified for multiple biogenic species
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk,rgasuniv   !constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal
        use canopy_phot_mod
        use canopy_bioparm_mod

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: ZK(:)           ! Input model heights (m)
        REAL(RK),    INTENT( IN )       :: FCLAI(:)        ! Input Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )       :: FCH             ! Model input canopy height (m)
        REAL(RK),    INTENT( IN )       :: LAI             ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )       :: CLU             ! Model input Clumping Index
        REAL(RK),    INTENT( IN )       :: COSZEN          ! Model input Cosine Solar Zenith Angle
        REAL(RK),    INTENT( IN )       :: SFCRAD          ! Model input Instantaneous surface downward shortwave flux (W/m2)
        REAL(RK),    INTENT( IN )       :: TEMP2           ! Model input 2-m Temperature (K)
        INTEGER,     INTENT( IN )       :: LU_OPT          ! integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
        INTEGER,     INTENT( IN )       :: VTYPE           ! Grid cell dominant vegetation type
        INTEGER,     INTENT( IN )       :: EMI_IND         ! Input biogenic emissions index
        REAL(RK),    INTENT( OUT )      :: EMI_OUT(:)      ! Output emissions (kg /m2 s)

! Local Variables
        REAL(RK) :: RJCF(SIZE(ZK))                 ! Photolysis correction factor for sun/shade fraction of leaf layer
        REAL(RK) :: PPFD_SUN(SIZE(ZK))             ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK) :: PPFD_SHADE(SIZE(ZK))           ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK) :: ATEMP_SUN(SIZE(ZK))            ! Regression coefficient A for sun leaves (Silva et al., 2020)
        REAL(RK) :: BTEMP_SUN(SIZE(ZK))            ! Regression coefficient B for sun leaves
        REAL(RK) :: ATEMP_SHADE(SIZE(ZK))          ! Regression coefficient A for shade leaves
        REAL(RK) :: BTEMP_SHADE(SIZE(ZK))          ! Regression coefficient B for shade leaves
        REAL(RK) :: TLEAF_SUN(SIZE(ZK))            ! Leaf temp for sunlit leaves (K)
        REAL(RK) :: TLEAF_SHADE(SIZE(ZK))          ! Leaf temp for shaded leaves (K)
        REAL(RK) :: TLEAF_AVE(SIZE(ZK))            ! Average Leaf temp (K)
        REAL(RK) :: TLEAF24_AVE(SIZE(ZK))          ! Average Leaf temp over the past 24 hours (K)
        REAL(RK) :: TLEAF240_AVE(SIZE(ZK))         ! Average Leaf temp over the past 240 hours (K)
        REAL(RK) :: GammaTLEAF_SUN_NUM(SIZE(ZK))   ! Numerator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE_NUM(SIZE(ZK)) ! Numerator in Tleaf shade activity factor
        REAL(RK) :: GammaTLEAF_SUN_DEN(SIZE(ZK))   ! Denominator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE_DEN(SIZE(ZK)) ! Denominator in Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SUN(SIZE(ZK))       ! Tleaf sun activity factor
        REAL(RK) :: GammaTLEAF_SHADE(SIZE(ZK))     ! Tleaf shade activity factor
        REAL(RK) :: GammaTLEAF_AVE(SIZE(ZK))       ! Average Tleaf activity factor
        REAL(RK) :: PPFD24_SUN(SIZE(ZK))           ! Average PPFD sun over the past 24 hours
        REAL(RK) :: PPFD24_SHADE(SIZE(ZK))         ! Average PPFD shade over the past 24 hours
        REAL(RK) :: PPFD240_SUN(SIZE(ZK))          ! Average PPFD sun over the past 240 hours
        REAL(RK) :: PPFD240_SHADE(SIZE(ZK))        ! Average PPFD shade over the past 240 hours
        REAL(RK) :: CP_SUN(SIZE(ZK))               ! Normalized emission capacity sun at PPFD = 1000 umol phot/m2 s
        REAL(RK) :: CP_SHADE(SIZE(ZK))             ! Normalized emission capacity shade at PPFD = 1000 umol phot/m2 s
        REAL(RK) :: ALPHA_P_SUN(SIZE(ZK))          ! Quantum yield of isoprene sunlit (mol/mol)
        REAL(RK) :: ALPHA_P_SHADE(SIZE(ZK))        ! Quantum yield of isoprene shade (mol/mol)
        REAL(RK) :: GammaPPFD_SUN(SIZE(ZK))        ! PPFD activity factor sun (unitless)
        REAL(RK) :: GammaPPFD_SHADE(SIZE(ZK))      ! PPFD activity factor shade (unitless)
        REAL(RK) :: GammaPPFD_AVE(SIZE(ZK))        ! PPFD activity factor ave sun and shade
        REAL(RK) :: E_OPT(SIZE(ZK))                ! maximum normalized emission capacity
        REAL(RK) :: TLEAF_OPT(SIZE(ZK))            ! Tleaf at which E_OPT occurs (K)
        REAL(RK) :: FLAI(SIZE(ZK))                 ! Fractional LAI in layer
        REAL(RK) :: CT1                            ! Activation energy (kJ/mol)
        REAL(RK) :: CEO                            ! Empirical coefficient
        REAL(RK) :: EF                             ! Final Mapped Emission factor (EF) (ug/m2 hr)
        integer i

! Constant Canopy Parameters
        REAL(RK),          PARAMETER     :: PPFD_CONST      =  4.5_rk     !based on MEGANv3 4.5 (umol photons/J)
        REAL(RK),          PARAMETER     :: PPFD0_SUN       =  200.0      !Constant PPFDo sunlit (umol/m2 s) (Guenther et al.,2012)
        REAL(RK),          PARAMETER     :: PPFD0_SHADE     =  50.0       !Constant PPFDo shaded (umol/m2 s) (Guenther et al.,2012)
        REAL(RK),          PARAMETER     :: ATEMP_TOP_SUN   =  -13.891_rk !Linearized 2-m temp --> leaf temp parameters
        REAL(RK),          PARAMETER     :: ATEMP_MID_SUN   =  -1.032_rk  !Based on Table 1 in Silva et al. (2022)
        REAL(RK),          PARAMETER     :: ATEMP_BOT_SUN   =  -5.589_rk  !...
        REAL(RK),          PARAMETER     :: BTEMP_TOP_SUN   =  1.064_rk   !...
        REAL(RK),          PARAMETER     :: BTEMP_MID_SUN   =  1.031_rk   !...
        REAL(RK),          PARAMETER     :: BTEMP_BOT_SUN   =  1.051_rk   !...
        REAL(RK),          PARAMETER     :: ATEMP_TOP_SHADE =  -12.846_rk !...
        REAL(RK),          PARAMETER     :: ATEMP_MID_SHADE =  -1.068_rk  !...
        REAL(RK),          PARAMETER     :: ATEMP_BOT_SHADE =  -5.955_rk  !...
        REAL(RK),          PARAMETER     :: BTEMP_TOP_SHADE =  1.060_rk   !...
        REAL(RK),          PARAMETER     :: BTEMP_MID_SHADE =  1.031_rk   !...
        REAL(RK),          PARAMETER     :: BTEMP_BOT_SHADE =  1.053_rk   !...
        REAL(RK),          PARAMETER     :: CT2             =  230.0_rk   !Deactivation energy (kJ/mol) (Guenther et al., 2012)


!Calculate photolyis shading/correction factor through canopy, i.e., the fraction of sunlit leaves downward through canopy

        call canopy_phot(FCLAI, LAI, CLU, COSZEN, RJCF)

! Calculate PPFD sun and shade through canopy layers using photolysis reduction factorsi (above canopy RJCF = 1)

        PPFD_SUN   = RJCF*SFCRAD*PPFD_CONST
        PPFD_SHADE = MAX((1.0-RJCF),0.0_rk)*SFCRAD*PPFD_CONST
!
! Use linear canopy temperature model based on Silva et al. (2020) to get approx. sun/shade leaf temperatures
! through canopy (ignores effect of wind speed on leaf boundary layer ~ 1 % error/bias)
!Citation:
!Silva, S. J., Heald, C. L., and Guenther, A. B.: Development of a reduced-complexity plant canopy
!physics surrogate model for use in chemical transport models: a case study with GEOS-Chem v12.3.0,
!Geosci. Model Dev., 13, 2569–2585, https://doi.org/10.5194/gmd-13-2569-2020, 2020.
        do i=1, SIZE(ZK)  !calculate linear change in parameters in top to mid and mid to bottom of canopy regions
            if (ZK(i) .gt. FCH) then ! above canopy, Tleaf = Tair
                ATEMP_SUN(i)   = 0.0
                BTEMP_SUN(i)   = 1.0
                ATEMP_SHADE(i) = 0.0
                BTEMP_SHADE(i) = 1.0
            else if (ZK(i) .le. FCH .and. ZK(i) .gt. FCH/2) then  !top to mid
                ATEMP_SUN(i)   = interp_linear1_internal((/ FCH/2,FCH /), (/ ATEMP_MID_SUN,ATEMP_TOP_SUN /),ZK(i))
                BTEMP_SUN(i)   = interp_linear1_internal((/ FCH/2,FCH /), (/ BTEMP_MID_SUN,BTEMP_TOP_SUN /),ZK(i))
                ATEMP_SHADE(i) = interp_linear1_internal((/ FCH/2,FCH /), (/ ATEMP_MID_SHADE,ATEMP_TOP_SHADE /),ZK(i))
                BTEMP_SHADE(i) = interp_linear1_internal((/ FCH/2,FCH /), (/ BTEMP_MID_SHADE,BTEMP_TOP_SHADE /),ZK(i))
            else if (ZK(i) .le. FCH/2) then !mid to bottom
                ATEMP_SUN(i)   = interp_linear1_internal((/ ZK(1),FCH/2 /), (/ ATEMP_BOT_SUN,ATEMP_MID_SUN /),ZK(i))
                BTEMP_SUN(i)   = interp_linear1_internal((/ ZK(1),FCH/2 /), (/ BTEMP_BOT_SUN,BTEMP_MID_SUN /),ZK(i))
                ATEMP_SHADE(i) = interp_linear1_internal((/ ZK(1),FCH/2 /), (/ ATEMP_BOT_SHADE,ATEMP_MID_SHADE /),ZK(i))
                BTEMP_SHADE(i) = interp_linear1_internal((/ ZK(1),FCH/2 /), (/ BTEMP_BOT_SHADE,BTEMP_MID_SHADE /),ZK(i))
            end if
        end do

        TLEAF_SUN   = ATEMP_SUN + (BTEMP_SUN*TEMP2)
        TLEAF_SHADE = ATEMP_SHADE + (BTEMP_SHADE*TEMP2)
        TLEAF_AVE = (TLEAF_SUN*RJCF) + (TLEAF_SHADE*(1.0-RJCF)) ! average = sum sun and shade weighted by sunlit fraction

! Calculate maximum normalized emission capacity (E_OPT) and Tleaf at E_OPT
        TLEAF240_AVE   = TLEAF_AVE  !Assume instantaneous TLEAF estimate to get TLEAF240 and TLEAF24 (improve...)
        TLEAF24_AVE    = TLEAF_AVE
        TLEAF_OPT = 313.0_rk + (0.6_rk * (TLEAF240_AVE-297.0_rk)) !Guenther et al. (2012)

! Calculate emission species/plant-dependent mapped emission factors
        call canopy_biop(EMI_IND, LU_OPT, VTYPE, EF, CT1, CEO)

        E_OPT = CEO * EXP(0.05_rk * (TLEAF24_AVE-297.0_rk)) * EXP(0.05_rk * (TLEAF240_AVE-297.0_rk))

! Calculate gamma (activity) values for average Tleaf (Clifton et al., 2022)
        GammaTLEAF_SUN_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN)))
        GammaTLEAF_SUN_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SUN))))
        GammaTLEAF_SUN     = E_OPT*(GammaTLEAF_SUN_NUM/GammaTLEAF_SUN_DEN)

        GammaTLEAF_SHADE_NUM = CT2*exp((CT1/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE)))
        GammaTLEAF_SHADE_DEN = (CT2-CT1)*(1.0-exp((CT2/(rgasuniv/1000.0))*((1.0/TLEAF_OPT)-(1.0/TLEAF_SHADE))))
        GammaTLEAF_SHADE     = E_OPT*(GammaTLEAF_SHADE_NUM/GammaTLEAF_SHADE_DEN)

        GammaTLEAF_AVE = (GammaTLEAF_SUN*RJCF) + (GammaTLEAF_SHADE*(1.0-RJCF)) ! average = sum sun and shade weighted by sunlit fraction

! Calculate gamma (activity) values for average PPFD (Clifton et al., 2022)
        PPFD240_SUN   = PPFD_SUN/2.0  !Clifton et al...halve the instantaneous PPFD estimate to get PPFD240 and PPFD24 (improve...)
        PPFD240_SHADE = PPFD_SHADE/2.0
        PPFD24_SUN    = PPFD_SUN/2.0
        PPFD24_SHADE  = PPFD_SHADE/2.0

        CP_SUN = 0.0468*(PPFD240_SUN**(0.6))*exp(0.005*(PPFD24_SUN-PPFD0_SUN))
        CP_SHADE = 0.0468*(PPFD240_SUN**(0.6))*exp(0.005*(PPFD24_SHADE-PPFD0_SHADE))
        ALPHA_P_SUN = 0.004 - 0.0005*log(PPFD240_SUN)
        ALPHA_P_SHADE = 0.004 - 0.0005*log(PPFD240_SHADE)
        GammaPPFD_SUN   = CP_SUN*((ALPHA_P_SUN*PPFD_SUN)/SQRT(1.0 + (ALPHA_P_SUN**2.0) * (PPFD_SUN**2.0)))
        GammaPPFD_SHADE = CP_SHADE*((ALPHA_P_SHADE*PPFD_SHADE)/SQRT(1.0 + (ALPHA_P_SHADE**2.0) * (PPFD_SHADE**2.0)))

        GammaPPFD_AVE = (GammaPPFD_SUN*RJCF) + (GammaPPFD_SHADE*(1.0-RJCF)) ! average = sum sun and shade weighted by sunlit fraction

! Calculate isoprene emissions profile in the canopy
        EMI_OUT = 0.0_rk  ! set initial emissions profile to zero
        do i=1, SIZE(ZK)
            if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then  ! above ground level and at/below canopy top
                FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/FCH    !fractional LAI scaled to canopy depth (~ LAD)
                EMI_OUT(i) = FLAI(i) * EF * GammaTLEAF_AVE(i) * GammaPPFD_AVE(i)  ! (ug m-2 hr-1)
                EMI_OUT(i) = EMI_OUT(i) * 2.7777777777778E-13_rk !TBD:  convert emissions output to (kg m-2 s-1)
            end if
        end do

    END SUBROUTINE CANOPY_BIO

end module canopy_bioemi_mod
