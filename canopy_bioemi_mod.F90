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
        REAL(RK) :: EF1,EF2,EF3,EF4,EF5,EF6,EF7    ! Plant Emission factors (EF) (ug/m2 hr)
        REAL(RK) :: EF8,EF9,EF10,EF11,EF12,EF13    ! Plant Emission factors (EF) (ug/m2 hr)
        REAL(RK) :: EF14,EF15                      ! Plant Emission factors (EF) (ug/m2 hr)
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

! Plant-Dependent emissions capacity/factors (EFs) for Isoprene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_ISOP    =  600.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_ISOP    =  3000.0_rk    ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_ISOP    =  1.0_rk       ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_ISOP    =  7000.0_rk    ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_ISOP    =  10000.0_rk   ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_ISOP    =  7000.0_rk    ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_ISOP    =  10000.0_rk   ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_ISOP    =  11000.0_rk   ! Broadleaf Deciduous Boreal Tree                                                                          
        REAL(RK),          PARAMETER     :: EF9_ISOP    =  2000.0_rk    ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_ISOP   =  4000.0_rk    ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_ISOP   =  4000.0_rk    ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_ISOP   =  1600.0_rk    ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_ISOP   =  800.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_ISOP   =  200.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_ISOP   =  1.0_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Myrcene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_MYRC    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_MYRC    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_MYRC    =  60.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_MYRC    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_MYRC    =  30.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_MYRC    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_MYRC    =  30.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_MYRC    =  30.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_MYRC    =  30.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_MYRC   =  50.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_MYRC   =  30.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_MYRC   =  0.3_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_MYRC   =  0.3_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_MYRC   =  0.3_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_MYRC   =  0.3_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Sabinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_SABI    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_SABI    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_SABI    =  40.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_SABI    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_SABI    =  50.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_SABI    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_SABI    =  50.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_SABI    =  50.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_SABI    =  50.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_SABI   =  70.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_SABI   =  50.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_SABI   =  0.7_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_SABI   =  0.7_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_SABI   =  0.7_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_SABI   =  0.7_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Limonene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_LIMO    =  100.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_LIMO    =  100.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_LIMO    =  130.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_LIMO    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_LIMO    =  80.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_LIMO    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_LIMO    =  80.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_LIMO    =  80.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_LIMO    =  60.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_LIMO   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_LIMO   =  60.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_LIMO   =  0.7_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_LIMO   =  0.7_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_LIMO   =  0.7_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_LIMO   =  0.7_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for 3-Carene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_CARE    =  160.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_CARE    =  160.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_CARE    =  80.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_CARE    =  40.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_CARE    =  30.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_CARE    =  40.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_CARE    =  30.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_CARE    =  30.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_CARE    =  30.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_CARE   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_CARE   =  30.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_CARE   =  0.3_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_CARE   =  0.3_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_CARE   =  0.3_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_CARE   =  0.3_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for t-Beta-Ocimene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_OCIM    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_OCIM    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_OCIM    =  60.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_OCIM    =  150.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_OCIM    =  120.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_OCIM    =  150.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_OCIM    =  120.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_OCIM    =  120.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_OCIM    =  90.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_OCIM   =  150.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_OCIM   =  90.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_OCIM   =  2.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_OCIM   =  2.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_OCIM   =  2.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_OCIM   =  2.0_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Beta-Pinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_BPIN    =  300.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_BPIN    =  300.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_BPIN    =  200.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_BPIN    =  120.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_BPIN    =  130.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_BPIN    =  120.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_BPIN    =  130.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_BPIN    =  130.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_BPIN    =  100.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_BPIN   =  150.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_BPIN   =  100.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_BPIN   =  1.5_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_BPIN   =  1.5_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_BPIN   =  1.5_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_BPIN   =  1.5_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Alpha-Pinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_APIN    =  500.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_APIN    =  500.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_APIN    =  510.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_APIN    =  600.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_APIN    =  400.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_APIN    =  600.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_APIN    =  400.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_APIN    =  400.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_APIN    =  200.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_APIN   =  300.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_APIN   =  200.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_APIN   =  2.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_APIN   =  2.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_APIN   =  2.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_APIN   =  2.0_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Other Monoterpenes (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! ! Other Monoterpenes category (34 compounds):  See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_MONO    =  180.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_MONO    =  180.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_MONO    =  170.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_MONO    =  150.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_MONO    =  150.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_MONO    =  150.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_MONO    =  150.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_MONO    =  150.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_MONO    =  110.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_MONO   =  200.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_MONO   =  110.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_MONO   =  5.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_MONO   =  5.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_MONO   =  5.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_MONO   =  5.0_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Alpha-Farnesene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_FARN    =  40.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_FARN    =  40.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_FARN    =  40.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_FARN    =  60.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_FARN    =  40.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_FARN    =  60.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_FARN    =  40.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_FARN    =  40.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_FARN    =  40.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_FARN   =  40.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_FARN   =  40.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_FARN   =  3.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_FARN   =  3.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_FARN   =  3.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_FARN   =  4.0_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Beta-Caryophyllene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_CARY    =  80.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_CARY    =  80.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_CARY    =  80.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_CARY    =  60.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_CARY    =  40.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_CARY    =  60.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_CARY    =  40.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_CARY    =  40.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_CARY    =  50.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_CARY   =  50.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_CARY   =  50.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_CARY   =  1.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_CARY   =  1.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_CARY   =  1.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_CARY   =  4.0_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Other Sesquieterpenes (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Other Sesquiterpenes category (30 compounds):  See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_SESQ    =  120.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_SESQ    =  120.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_SESQ    =  120.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_SESQ    =  120.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_SESQ    =  100.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_SESQ    =  120.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_SESQ    =  100.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_SESQ    =  100.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_SESQ    =  100.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_SESQ   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_SESQ   =  100.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_SESQ   =  1.0_rk       ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_SESQ   =  1.0_rk       ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_SESQ   =  1.0_rk       ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_SESQ   =  1.0_rk       ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for 232-MBO (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_MBOL    =  700.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_MBOL    =  60.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_MBOL    =  0.01_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_MBOL    =  0.01_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_MBOL    =  0.01_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_MBOL    =  0.01_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_MBOL    =  0.01_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_MBOL    =  2.0_rk       ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_MBOL    =  0.01_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_MBOL   =  0.01_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_MBOL   =  0.01_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_MBOL   =  0.01_rk      ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_MBOL   =  0.01_rk      ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_MBOL   =  0.01_rk      ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_MBOL   =  0.01_rk      ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Methanol (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_METH    =  900.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_METH    =  900.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_METH    =  900.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_METH    =  500.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_METH    =  900.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_METH    =  500.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_METH    =  900.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_METH    =  900.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_METH    =  900.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_METH   =  900.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_METH   =  900.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_METH   =  500.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_METH   =  500.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_METH   =  500.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_METH   =  900.0_rk     ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Acetone (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_ACET    =  240.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_ACET    =  240.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_ACET    =  240.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_ACET    =  240.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_ACET    =  240.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_ACET    =  240.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_ACET    =  240.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_ACET    =  240.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_ACET    =  240.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_ACET   =  240.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_ACET   =  240.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_ACET   =  80.0_rk      ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_ACET   =  80.0_rk      ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_ACET   =  80.0_rk      ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_ACET   =  80.0_rk      ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Carbon Monoxide (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1_CO      =  600.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_CO      =  600.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_CO      =  600.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_CO      =  600.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_CO      =  600.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_CO      =  600.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_CO      =  600.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_CO      =  600.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_CO      =  600.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_CO     =  600.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_CO     =  600.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_CO     =  600.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_CO     =  600.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_CO     =  600.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_CO     =  600.0_rk     ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for BIDI VOC species (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Bidirectional VOC (5 compounds): See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_BVOC    =  500.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_BVOC    =  500.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_BVOC    =  500.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_BVOC    =  500.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_BVOC    =  500.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_BVOC    =  500.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_BVOC    =  500.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_BVOC    =  500.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_BVOC    =  500.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_BVOC   =  500.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_BVOC   =  500.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_BVOC   =  80.0_rk      ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_BVOC   =  80.0_rk      ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_BVOC   =  80.0_rk      ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_BVOC   =  80.0_rk      ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Stress VOCs (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Stress VOC (15 compounds): See Table 1 of Guenther et al. (2012) 
        REAL(RK),          PARAMETER     :: EF1_SVOC    =  300.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_SVOC    =  300.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_SVOC    =  300.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_SVOC    =  300.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_SVOC    =  300.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_SVOC    =  300.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_SVOC    =  300.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_SVOC    =  300.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_SVOC    =  300.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_SVOC   =  300.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_SVOC   =  300.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_SVOC   =  300.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_SVOC   =  300.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_SVOC   =  300.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_SVOC   =  300.0_rk     ! Crop1
! Plant-Dependent emissions capacity/factors (EFs) for Other VOCs (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Other VOC (49 compounds): See Table 1 of Guenther et al. (2012)
        REAL(RK),          PARAMETER     :: EF1_OVOC    =  140.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2_OVOC    =  140.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3_OVOC    =  140.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4_OVOC    =  140.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5_OVOC    =  140.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6_OVOC    =  140.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7_OVOC    =  140.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8_OVOC    =  140.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF9_OVOC    =  140.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10_OVOC   =  140.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11_OVOC   =  140.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12_OVOC   =  140.0_rk     ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13_OVOC   =  140.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14_OVOC   =  140.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15_OVOC   =  140.0_rk     ! Crop1
        
! Species-Dependent Parameterized Canopy Model Parameters (Table 4 of Guenther et al., 2012)
        REAL(RK),          PARAMETER     :: CT1_ISOP         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_ISOP         =  2.0_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_MYRC         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_MYRC         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_SABI         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_SABI         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_LIMO         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_LIMO         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_CARE         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_CARE         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_OCIM         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_OCIM         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_BPIN         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_BPIN         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_APIN         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_APIN         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_MONO         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_MONO         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_FARN         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_FARN         =  2.37_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_CARY         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_CARY         =  2.37_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_SESQ         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_SESQ         =  2.37_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_MBOL         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_MBOL         =  2.0_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_METH         =  60.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_METH         =  1.6_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_ACET         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_ACET         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_CO           =  60.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_CO           =  1.6_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_BVOC         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_BVOC         =  2.0_rk     !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_SVOC         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_SVOC         =  1.83_rk    !Empirical coefficient
        REAL(RK),          PARAMETER     :: CT1_OVOC         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_OVOC         =  1.83_rk    !Empirical coefficient


!Calculate photolyis shading/correction factor through canopy, i.e., the fraction of sunlit leaves downward through canopy

       call canopy_phot(FCLAI, LAI, CLU, COSZEN, RJCF)

        ! Calculate PPFD sun and shade through canopy layers using photolysis reduction factors
        ! above canopy RJCF = 1 

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

! Set tree and species dependent coefficients
       if (EMI_IND .eq. 1 ) then
          CT1  = CT1_ISOP
          CEO  = CEO_ISOP
          EF1  = EF1_ISOP
          EF2  = EF2_ISOP
          EF3  = EF3_ISOP
          EF4  = EF4_ISOP
          EF5  = EF5_ISOP
          EF6  = EF6_ISOP
          EF7  = EF7_ISOP
          EF8  = EF8_ISOP
          EF9  = EF9_ISOP
          EF10 = EF10_ISOP
          EF11 = EF11_ISOP
          EF12 = EF12_ISOP
          EF13 = EF13_ISOP
          EF14 = EF14_ISOP
          EF15 = EF15_ISOP 
       else if (EMI_IND .eq. 2 ) then 
          CT1 = CT1_MYRC
          CEO = CEO_MYRC
          EF1  = EF1_MYRC
          EF2  = EF2_MYRC
          EF3  = EF3_MYRC
          EF4  = EF4_MYRC
          EF5  = EF5_MYRC
          EF6  = EF6_MYRC
          EF7  = EF7_MYRC
          EF8  = EF8_MYRC
          EF9  = EF9_MYRC
          EF10 = EF10_MYRC
          EF11 = EF11_MYRC
          EF12 = EF12_MYRC
          EF13 = EF13_MYRC
          EF14 = EF14_MYRC
          EF15 = EF15_MYRC
       else if (EMI_IND .eq. 3 ) then 
          CT1 = CT1_SABI
          CEO = CEO_SABI
          EF1  = EF1_SABI
          EF2  = EF2_SABI
          EF3  = EF3_SABI
          EF4  = EF4_SABI
          EF5  = EF5_SABI
          EF6  = EF6_SABI
          EF7  = EF7_SABI
          EF8  = EF8_SABI
          EF9  = EF9_SABI
          EF10 = EF10_SABI
          EF11 = EF11_SABI
          EF12 = EF12_SABI
          EF13 = EF13_SABI
          EF14 = EF14_SABI
          EF15 = EF15_SABI
       else if (EMI_IND .eq. 4 ) then
          CT1 = CT1_LIMO
          CEO = CEO_LIMO
          EF1  = EF1_LIMO
          EF2  = EF2_LIMO
          EF3  = EF3_LIMO
          EF4  = EF4_LIMO
          EF5  = EF5_LIMO
          EF6  = EF6_LIMO
          EF7  = EF7_LIMO
          EF8  = EF8_LIMO
          EF9  = EF9_LIMO
          EF10 = EF10_LIMO
          EF11 = EF11_LIMO
          EF12 = EF12_LIMO
          EF13 = EF13_LIMO
          EF14 = EF14_LIMO
          EF15 = EF15_LIMO
       else if (EMI_IND .eq. 5 ) then
          CT1 = CT1_CARE
          CEO = CEO_CARE
          EF1  = EF1_CARE
          EF2  = EF2_CARE
          EF3  = EF3_CARE
          EF4  = EF4_CARE
          EF5  = EF5_CARE
          EF6  = EF6_CARE
          EF7  = EF7_CARE
          EF8  = EF8_CARE
          EF9  = EF9_CARE
          EF10 = EF10_CARE
          EF11 = EF11_CARE
          EF12 = EF12_CARE
          EF13 = EF13_CARE
          EF14 = EF14_CARE
          EF15 = EF15_CARE
       else if (EMI_IND .eq. 6 ) then
          CT1 = CT1_OCIM
          CEO = CEO_OCIM
          EF1  = EF1_OCIM
          EF2  = EF2_OCIM
          EF3  = EF3_OCIM
          EF4  = EF4_OCIM
          EF5  = EF5_OCIM
          EF6  = EF6_OCIM
          EF7  = EF7_OCIM
          EF8  = EF8_OCIM
          EF9  = EF9_OCIM
          EF10 = EF10_OCIM
          EF11 = EF11_OCIM
          EF12 = EF12_OCIM
          EF13 = EF13_OCIM
          EF14 = EF14_OCIM
          EF15 = EF15_OCIM
       else if (EMI_IND .eq. 7 ) then
          CT1 = CT1_BPIN
          CEO = CEO_BPIN
          EF1  = EF1_BPIN
          EF2  = EF2_BPIN
          EF3  = EF3_BPIN
          EF4  = EF4_BPIN
          EF5  = EF5_BPIN
          EF6  = EF6_BPIN
          EF7  = EF7_BPIN
          EF8  = EF8_BPIN
          EF9  = EF9_BPIN
          EF10 = EF10_BPIN
          EF11 = EF11_BPIN
          EF12 = EF12_BPIN
          EF13 = EF13_BPIN
          EF14 = EF14_BPIN
          EF15 = EF15_BPIN
       else if (EMI_IND .eq. 8 ) then
          CT1 = CT1_APIN
          CEO = CEO_APIN
          EF1  = EF1_APIN
          EF2  = EF2_APIN
          EF3  = EF3_APIN
          EF4  = EF4_APIN
          EF5  = EF5_APIN
          EF6  = EF6_APIN
          EF7  = EF7_APIN
          EF8  = EF8_APIN
          EF9  = EF9_APIN
          EF10 = EF10_APIN
          EF11 = EF11_APIN
          EF12 = EF12_APIN
          EF13 = EF13_APIN
          EF14 = EF14_APIN
          EF15 = EF15_APIN
       else if (EMI_IND .eq. 9 ) then
          CT1 = CT1_MONO
          CEO = CEO_MONO
          EF1  = EF1_MONO
          EF2  = EF2_MONO
          EF3  = EF3_MONO
          EF4  = EF4_MONO
          EF5  = EF5_MONO
          EF6  = EF6_MONO
          EF7  = EF7_MONO
          EF8  = EF8_MONO
          EF9  = EF9_MONO
          EF10 = EF10_MONO
          EF11 = EF11_MONO
          EF12 = EF12_MONO
          EF13 = EF13_MONO
          EF14 = EF14_MONO
          EF15 = EF15_MONO
       else if (EMI_IND .eq. 10 ) then
          CT1 = CT1_FARN
          CEO = CEO_FARN
          EF1  = EF1_FARN
          EF2  = EF2_FARN
          EF3  = EF3_FARN
          EF4  = EF4_FARN
          EF5  = EF5_FARN
          EF6  = EF6_FARN
          EF7  = EF7_FARN
          EF8  = EF8_FARN
          EF9  = EF9_FARN
          EF10 = EF10_FARN
          EF11 = EF11_FARN
          EF12 = EF12_FARN
          EF13 = EF13_FARN
          EF14 = EF14_FARN
          EF15 = EF15_FARN
       else if (EMI_IND .eq. 11 ) then
          CT1 = CT1_CARY
          CEO = CEO_CARY
          EF1  = EF1_CARY
          EF2  = EF2_CARY
          EF3  = EF3_CARY
          EF4  = EF4_CARY
          EF5  = EF5_CARY
          EF6  = EF6_CARY
          EF7  = EF7_CARY
          EF8  = EF8_CARY
          EF9  = EF9_CARY
          EF10 = EF10_CARY
          EF11 = EF11_CARY
          EF12 = EF12_CARY
          EF13 = EF13_CARY
          EF14 = EF14_CARY
          EF15 = EF15_CARY
       else if (EMI_IND .eq. 12 ) then
          CT1 = CT1_SESQ
          CEO = CEO_SESQ
          EF1  = EF1_SESQ
          EF2  = EF2_SESQ
          EF3  = EF3_SESQ
          EF4  = EF4_SESQ
          EF5  = EF5_SESQ
          EF6  = EF6_SESQ
          EF7  = EF7_SESQ
          EF8  = EF8_SESQ
          EF9  = EF9_SESQ
          EF10 = EF10_SESQ
          EF11 = EF11_SESQ
          EF12 = EF12_SESQ
          EF13 = EF13_SESQ
          EF14 = EF14_SESQ
          EF15 = EF15_SESQ
       else if (EMI_IND .eq. 13 ) then
          CT1 = CT1_MBOL
          CEO = CEO_MBOL
          EF1  = EF1_MBOL
          EF2  = EF2_MBOL
          EF3  = EF3_MBOL
          EF4  = EF4_MBOL
          EF5  = EF5_MBOL
          EF6  = EF6_MBOL
          EF7  = EF7_MBOL
          EF8  = EF8_MBOL
          EF9  = EF9_MBOL
          EF10 = EF10_MBOL
          EF11 = EF11_MBOL
          EF12 = EF12_MBOL
          EF13 = EF13_MBOL
          EF14 = EF14_MBOL
          EF15 = EF15_MBOL
       else if (EMI_IND .eq. 14 ) then
          CT1 = CT1_METH
          CEO = CEO_METH
          EF1  = EF1_METH
          EF2  = EF2_METH
          EF3  = EF3_METH
          EF4  = EF4_METH
          EF5  = EF5_METH
          EF6  = EF6_METH
          EF7  = EF7_METH
          EF8  = EF8_METH
          EF9  = EF9_METH
          EF10 = EF10_METH
          EF11 = EF11_METH
          EF12 = EF12_METH
          EF13 = EF13_METH
          EF14 = EF14_METH
          EF15 = EF15_METH
       else if (EMI_IND .eq. 15 ) then
          CT1 = CT1_ACET
          CEO = CEO_ACET
          EF1  = EF1_ACET
          EF2  = EF2_ACET
          EF3  = EF3_ACET
          EF4  = EF4_ACET
          EF5  = EF5_ACET
          EF6  = EF6_ACET
          EF7  = EF7_ACET
          EF8  = EF8_ACET
          EF9  = EF9_ACET
          EF10 = EF10_ACET
          EF11 = EF11_ACET
          EF12 = EF12_ACET
          EF13 = EF13_ACET
          EF14 = EF14_ACET
          EF15 = EF15_ACET
       else if (EMI_IND .eq. 16 ) then
          CT1 = CT1_CO
          CEO = CEO_CO
          EF1  = EF1_CO
          EF2  = EF2_CO
          EF3  = EF3_CO
          EF4  = EF4_CO
          EF5  = EF5_CO
          EF6  = EF6_CO
          EF7  = EF7_CO
          EF8  = EF8_CO
          EF9  = EF9_CO
          EF10 = EF10_CO
          EF11 = EF11_CO
          EF12 = EF12_CO
          EF13 = EF13_CO
          EF14 = EF14_CO
          EF15 = EF15_CO
       else if (EMI_IND .eq. 17 ) then
          CT1 = CT1_BVOC
          CEO = CEO_BVOC
          EF1  = EF1_BVOC
          EF2  = EF2_BVOC
          EF3  = EF3_BVOC
          EF4  = EF4_BVOC
          EF5  = EF5_BVOC
          EF6  = EF6_BVOC
          EF7  = EF7_BVOC
          EF8  = EF8_BVOC
          EF9  = EF9_BVOC
          EF10 = EF10_BVOC
          EF11 = EF11_BVOC
          EF12 = EF12_BVOC
          EF13 = EF13_BVOC
          EF14 = EF14_BVOC
          EF15 = EF15_BVOC
       else if (EMI_IND .eq. 18 ) then
          CT1 = CT1_SVOC
          CEO = CEO_SVOC
          EF1  = EF1_SVOC
          EF2  = EF2_SVOC
          EF3  = EF3_SVOC
          EF4  = EF4_SVOC
          EF5  = EF5_SVOC
          EF6  = EF6_SVOC
          EF7  = EF7_SVOC
          EF8  = EF8_SVOC
          EF9  = EF9_SVOC
          EF10 = EF10_SVOC
          EF11 = EF11_SVOC
          EF12 = EF12_SVOC
          EF13 = EF13_SVOC
          EF14 = EF14_SVOC
          EF15 = EF15_SVOC
       else   ! EMI_IND = 19 
          CT1 = CT1_OVOC
          CEO = CEO_OVOC
          EF1  = EF1_OVOC
          EF2  = EF2_OVOC
          EF3  = EF3_OVOC
          EF4  = EF4_OVOC
          EF5  = EF5_OVOC
          EF6  = EF6_OVOC
          EF7  = EF7_OVOC
          EF8  = EF8_OVOC
          EF9  = EF9_OVOC
          EF10 = EF10_OVOC
          EF11 = EF11_OVOC
          EF12 = EF12_OVOC
          EF13 = EF13_OVOC
          EF14 = EF14_OVOC
          EF15 = EF15_OVOC
       end if

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

       if (LU_OPT .eq. 0 .or. LU_OPT .eq. 1) then !VIIRS or MODIS  LU types

! Simple MEGAN (Table 3 in Guenther et al., 2012) PFT to VIIRS/MODIS VTYPE mapping
       if (VTYPE .eq. 1) then !VIIRS Cat 1 Evergreen Needleleaf
        !--> Average Needleleaf Evergreen Temperate Tree and Needleleaf Evergreen Boreal Tree

            EF = (EF1+EF2)/2.0_rk

        else if (VTYPE .eq. 2) then !VIIRS/MODIS Cat 2 Evergreen Broadleaf
        !--> Average Broadleaf Evergreen Tropical Tree and Broadleaf Evergreen Temperate Tree

            EF = (EF4+EF5)/2.0_rk

        else if (VTYPE .eq. 3) then !VIIRS/MODIS Cat 3 Deciduous Needleaf
        !--> Average Needleleaf Deciduous Boreal Tree

            EF = EF3

        else if (VTYPE .eq. 4) then !VIIRS/MODIS Cat 4 Deciduous Broadleaf
        !--> Average Broadleaf Deciduous Tropical Tree, Broadleaf Deciduous Temperate Tree, 
        ! and Broadleaf Deciduous Boreal Tree

            EF = (EF6+EF7+EF8)/3.0_rk

        else if (VTYPE .eq. 5) then !VIIRS/MODIS Cat 5 Mixed forests
        !--> Avearge of all above EF1-EF8 PFTs.

            EF = (EF1+EF2+EF3+EF4+EF5+EF6+EF7+EF8)/8.0_rk

        else if (VTYPE .ge. 6 .and. VTYPE .le. 7) then !VIIRS/MODIS Cat 6-7 Closed/Open Shrublands
        !--> Avearge Broadleaf Evergreen Temperate Shrub, Broadleaf Deciduous Temperate Shrub,
        ! and Broadleaf Deciduous Boreal Shrub

            EF = (EF9+EF10+EF11)/3.0_rk

        else if (VTYPE .ge. 8 .and. VTYPE .le. 10) then !VIIRS/MODIS Cat 8-10 Savannas and Grasslands
        !--> Avearge Arctic C3 Grass, Cool C3 Grass, Warm C4 Grass)

            EF = (EF12+EF13+EF14)/3.0_rk

        else if (VTYPE .eq. 12 ) then !VIIRS/MODIS Cat 12 Croplands
        !--> Crop1

            EF = EF15

        else 

            EF = 0.0_rk  !set EF=0 for non-valid VTYPEs

        end if

        else
            write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS (0) or MODIS (1) available right now...exiting'
            call exit(2)
        end if

! TBD:  Add input vegtype from model and select representative plant functional type for EPS/EF factors...
! Calculate isoprene emissions profile in the canopy
       EMI_OUT = 0.0_rk  ! set initial emissions profile to zero
       do i=1, SIZE(ZK)
           if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then  ! above ground level and at/below canopy top
               FLAI(i) = ((FCLAI(i+1) - FCLAI(i)) * LAI)/FCH    !fractional LAI scaled to canopy depth (~ LAD)
               EMI_OUT(i) = FLAI(i) * EF * GammaTLEAF_AVE(i) * GammaPPFD_AVE(i)  ! (ug/m2 hr)
               EMI_OUT(i) = EMI_OUT(i) * 2.77778E-13 !TBD:  convert emissions output to kg/m2 s
           end if
       end do 

    END SUBROUTINE CANOPY_BIO

end module canopy_bioemi_mod
