

# File canopy\_bioparm\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_bioparm\_mod.F90**](canopy__bioparm__mod_8F90.md)

[Go to the documentation of this file](canopy__bioparm__mod_8F90.md)


```Fortran



module canopy_bioparm_mod

    implicit none

contains

    SUBROUTINE canopy_biop( EMI_IND, LU_OPT, VTYPE, &
        EF, LDF, BETA, CT1, CEO, ANEW, AGRO, AMAT, AOLD, &
        ROOTA, ROOTB, CAQ, TAQ, DTAQ, CHT, THT, DTHT, &
        CLT, TLT, DTLT, CHW, THW, DTHW)

        use canopy_const_mod, ONLY: rk

        INTEGER,     INTENT( IN )       :: emi_ind
        INTEGER,     INTENT( IN )       :: lu_opt
        INTEGER,     INTENT( IN )       :: vtype

        REAL(rk),    INTENT( OUT )      :: ef
        REAL(rk),    INTENT( OUT )      :: ldf
        REAL(rk),    INTENT( OUT )      :: beta
        REAL(rk),    INTENT( OUT )      :: ct1
        REAL(rk),    INTENT( OUT )      :: ceo
        REAL(rk),    INTENT( OUT )      :: anew, agro, amat, aold
        REAL(rk),    INTENT( OUT )      :: roota, rootb
        REAL(rk),    INTENT( OUT )      :: caq
        REAL(rk),    INTENT( OUT )      :: taq
        REAL(rk),    INTENT( OUT )      :: dtaq
        REAL(rk),    INTENT( OUT )      :: cht
        REAL(rk),    INTENT( OUT )      :: tht
        REAL(rk),    INTENT( OUT )      :: dtht
        REAL(rk),    INTENT( OUT )      :: clt
        REAL(rk),    INTENT( OUT )      :: tlt
        REAL(rk),    INTENT( OUT )      :: dtlt
        REAL(rk),    INTENT( OUT )      :: chw
        REAL(rk),    INTENT( OUT )      :: thw
        REAL(rk),    INTENT( OUT )      :: dthw

        REAL(rk) :: ef1,ef2,ef3,ef4,ef5,ef6,ef7
        REAL(rk) :: ef8,ef9,ef10,ef11,ef12,ef13
        REAL(rk) :: ef14,ef15


        REAL(rk),          PARAMETER     :: ef1_isop    =  600.0_rk
        REAL(rk),          PARAMETER     :: ef2_isop    =  3000.0_rk
        REAL(rk),          PARAMETER     :: ef3_isop    =  1.0_rk
        REAL(rk),          PARAMETER     :: ef4_isop    =  7000.0_rk
        REAL(rk),          PARAMETER     :: ef5_isop    =  10000.0_rk
        REAL(rk),          PARAMETER     :: ef6_isop    =  7000.0_rk
        REAL(rk),          PARAMETER     :: ef7_isop    =  10000.0_rk
        REAL(rk),          PARAMETER     :: ef8_isop    =  11000.0_rk
        REAL(rk),          PARAMETER     :: ef9_isop    =  2000.0_rk
        REAL(rk),          PARAMETER     :: ef10_isop   =  4000.0_rk
        REAL(rk),          PARAMETER     :: ef11_isop   =  4000.0_rk
        REAL(rk),          PARAMETER     :: ef12_isop   =  1600.0_rk
        REAL(rk),          PARAMETER     :: ef13_isop   =  800.0_rk
        REAL(rk),          PARAMETER     :: ef14_isop   =  200.0_rk
        REAL(rk),          PARAMETER     :: ef15_isop   =  1.0_rk

        REAL(rk),          PARAMETER     :: anew_isop  = 0.05_rk
        REAL(rk),          PARAMETER     :: agro_isop  = 0.6_rk
        REAL(rk),          PARAMETER     :: amat_isop  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_isop  = 0.9_rk


! Plant-Dependent emissions capacity/factors (EFs) for Myrcene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_myrc    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_myrc    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_myrc    =  60.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_myrc    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_myrc    =  30.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_myrc    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_myrc    =  30.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_myrc    =  30.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_myrc    =  30.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_myrc   =  50.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_myrc   =  30.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_myrc   =  0.3_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_myrc   =  0.3_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_myrc   =  0.3_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_myrc   =  0.3_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Myrcene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_myrc  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_myrc  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_myrc  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_myrc  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Sabinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_sabi    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_sabi    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_sabi    =  40.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_sabi    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_sabi    =  50.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_sabi    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_sabi    =  50.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_sabi    =  50.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_sabi    =  50.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_sabi   =  70.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_sabi   =  50.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_sabi   =  0.7_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_sabi   =  0.7_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_sabi   =  0.7_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_sabi   =  0.7_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Sabinene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_sabi  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_sabi  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_sabi  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_sabi  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Limonene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_limo    =  100.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_limo    =  100.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_limo    =  130.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_limo    =  80.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_limo    =  80.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_limo    =  80.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_limo    =  80.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_limo    =  80.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_limo    =  60.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_limo   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_limo   =  60.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_limo   =  0.7_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_limo   =  0.7_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_limo   =  0.7_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_limo   =  0.7_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Limonene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_limo  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_limo  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_limo  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_limo  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for 3-Carene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_care    =  160.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_care    =  160.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_care    =  80.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_care    =  40.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_care    =  30.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_care    =  40.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_care    =  30.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_care    =  30.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_care    =  30.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_care   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_care   =  30.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_care   =  0.3_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_care   =  0.3_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_care   =  0.3_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_care   =  0.3_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for 3-Carene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_care  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_care  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_care  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_care  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for t-Beta-Ocimene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_ocim    =  70.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_ocim    =  70.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_ocim    =  60.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_ocim    =  150.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_ocim    =  120.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_ocim    =  150.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_ocim    =  120.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_ocim    =  120.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_ocim    =  90.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_ocim   =  150.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_ocim   =  90.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_ocim   =  2.0_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_ocim   =  2.0_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_ocim   =  2.0_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_ocim   =  2.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for t-Beta-Ocimene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_ocim  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_ocim  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_ocim  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_ocim  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Beta-Pinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_bpin    =  300.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_bpin    =  300.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_bpin    =  200.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_bpin    =  120.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_bpin    =  130.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_bpin    =  120.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_bpin    =  130.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_bpin    =  130.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_bpin    =  100.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_bpin   =  150.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_bpin   =  100.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_bpin   =  1.5_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_bpin   =  1.5_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_bpin   =  1.5_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_bpin   =  1.5_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Beta-Pinene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_bpin  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_bpin  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_bpin  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_bpin  = 1.05_rk


! Plant-Dependent emissions capacity/factors (EFs) for Alpha-Pinene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_apin    =  500.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_apin    =  500.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_apin    =  510.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_apin    =  600.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_apin    =  400.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_apin    =  600.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_apin    =  400.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_apin    =  400.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_apin    =  200.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_apin   =  300.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_apin   =  200.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_apin   =  2.0_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_apin   =  2.0_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_apin   =  2.0_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_apin   =  2.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Alpha-Pinene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_apin  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_apin  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_apin  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_apin  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Other Monoterpenes (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! ! Other Monoterpenes category (34 compounds):  See Table 1 of Guenther et al. (2012)
        REAL(rk),          PARAMETER     :: ef1_mono    =  180.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_mono    =  180.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_mono    =  170.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_mono    =  150.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_mono    =  150.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_mono    =  150.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_mono    =  150.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_mono    =  150.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_mono    =  110.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_mono   =  200.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_mono   =  110.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_mono   =  5.0_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_mono   =  5.0_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_mono   =  5.0_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_mono   =  5.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Other Monoterpenes as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_mono  = 2.0_rk
        REAL(rk),          PARAMETER     :: agro_mono  = 1.8_rk
        REAL(rk),          PARAMETER     :: amat_mono  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_mono  = 1.05_rk

! Plant-Dependent emissions capacity/factors (EFs) for Alpha-Farnesene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_farn    =  40.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_farn    =  40.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_farn    =  40.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_farn    =  60.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_farn    =  40.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_farn    =  60.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_farn    =  40.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_farn    =  40.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_farn    =  40.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_farn   =  40.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_farn   =  40.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_farn   =  3.0_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_farn   =  3.0_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_farn   =  3.0_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_farn   =  4.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Alpha-Farnesene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_farn  = 0.4_rk
        REAL(rk),          PARAMETER     :: agro_farn  = 0.6_rk
        REAL(rk),          PARAMETER     :: amat_farn  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_farn  = 0.95_rk

! Plant-Dependent emissions capacity/factors (EFs) for Beta-Caryophyllene (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_cary    =  80.0_rk      ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_cary    =  80.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_cary    =  80.0_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_cary    =  60.0_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_cary    =  40.0_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_cary    =  60.0_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_cary    =  40.0_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_cary    =  40.0_rk      ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_cary    =  50.0_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_cary   =  50.0_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_cary   =  50.0_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_cary   =  1.0_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_cary   =  1.0_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_cary   =  1.0_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_cary   =  4.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Beta-Caryophyllene as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_cary  = 0.4_rk
        REAL(rk),          PARAMETER     :: agro_cary  = 0.6_rk
        REAL(rk),          PARAMETER     :: amat_cary  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_cary  = 0.95_rk

! Plant-Dependent emissions capacity/factors (EFs) for Other Sesquieterpenes (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Other Sesquiterpenes category (30 compounds):  See Table 1 of Guenther et al. (2012)
        REAL(rk),          PARAMETER     :: ef1_sesq    =  120.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_sesq    =  120.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_sesq    =  120.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_sesq    =  120.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_sesq    =  100.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_sesq    =  120.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_sesq    =  100.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_sesq    =  100.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_sesq    =  100.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_sesq   =  100.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_sesq   =  100.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_sesq   =  1.0_rk       ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_sesq   =  1.0_rk       ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_sesq   =  1.0_rk       ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_sesq   =  1.0_rk       ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Other Sesquieterpenes as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_sesq  = 0.4_rk
        REAL(rk),          PARAMETER     :: agro_sesq  = 0.6_rk
        REAL(rk),          PARAMETER     :: amat_sesq  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_sesq  = 0.95_rk

! Plant-Dependent emissions capacity/factors (EFs) for 232-MBO (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_mbol    =  700.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_mbol    =  60.0_rk      ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_mbol    =  0.01_rk      ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_mbol    =  0.01_rk      ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_mbol    =  0.01_rk      ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_mbol    =  0.01_rk      ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_mbol    =  0.01_rk      ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_mbol    =  2.0_rk       ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_mbol    =  0.01_rk      ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_mbol   =  0.01_rk      ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_mbol   =  0.01_rk      ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_mbol   =  0.01_rk      ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_mbol   =  0.01_rk      ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_mbol   =  0.01_rk      ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_mbol   =  0.01_rk      ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for 232-MBO as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_mbol  = 0.05_rk
        REAL(rk),          PARAMETER     :: agro_mbol  = 0.6_rk
        REAL(rk),          PARAMETER     :: amat_mbol  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_mbol  = 0.9_rk

! Plant-Dependent emissions capacity/factors (EFs) for Methanol (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_meth    =  900.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_meth    =  900.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_meth    =  900.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_meth    =  500.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_meth    =  900.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_meth    =  500.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_meth    =  900.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_meth    =  900.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_meth    =  900.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_meth   =  900.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_meth   =  900.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_meth   =  500.0_rk     ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_meth   =  500.0_rk     ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_meth   =  500.0_rk     ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_meth   =  900.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Methanol as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_meth  = 3.5_rk
        REAL(rk),          PARAMETER     :: agro_meth  = 3.0_rk
        REAL(rk),          PARAMETER     :: amat_meth  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_meth  = 1.2_rk

! Plant-Dependent emissions capacity/factors (EFs) for Acetone (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_acet    =  240.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_acet    =  240.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_acet    =  240.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_acet    =  240.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_acet    =  240.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_acet    =  240.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_acet    =  240.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_acet    =  240.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_acet    =  240.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_acet   =  240.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_acet   =  240.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_acet   =  80.0_rk      ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_acet   =  80.0_rk      ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_acet   =  80.0_rk      ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_acet   =  80.0_rk      ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Acetone as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_acet  = 1.0_rk
        REAL(rk),          PARAMETER     :: agro_acet  = 1.0_rk
        REAL(rk),          PARAMETER     :: amat_acet  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_acet  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for Carbon Monoxide (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(rk),          PARAMETER     :: ef1_co      =  600.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_co      =  600.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_co      =  600.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_co      =  600.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_co      =  600.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_co      =  600.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_co      =  600.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_co      =  600.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_co      =  600.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_co     =  600.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_co     =  600.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_co     =  600.0_rk     ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_co     =  600.0_rk     ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_co     =  600.0_rk     ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_co     =  600.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Carbon Monoxide as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_co  = 1.0_rk
        REAL(rk),          PARAMETER     :: agro_co  = 1.0_rk
        REAL(rk),          PARAMETER     :: amat_co  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_co  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for BIDI VOC species (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Bidirectional VOC (5 compounds): See Table 1 of Guenther et al. (2012)
        REAL(rk),          PARAMETER     :: ef1_bvoc    =  500.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_bvoc    =  500.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_bvoc    =  500.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_bvoc    =  500.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_bvoc    =  500.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_bvoc    =  500.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_bvoc    =  500.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_bvoc    =  500.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_bvoc    =  500.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_bvoc   =  500.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_bvoc   =  500.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_bvoc   =  80.0_rk      ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_bvoc   =  80.0_rk      ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_bvoc   =  80.0_rk      ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_bvoc   =  80.0_rk      ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for BIDI VOC species as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_bvoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: agro_bvoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: amat_bvoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_bvoc  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for Stress VOCs (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Stress VOC (15 compounds): See Table 1 of Guenther et al. (2012)
        REAL(rk),          PARAMETER     :: ef1_svoc    =  300.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_svoc    =  300.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_svoc    =  300.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_svoc    =  300.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_svoc    =  300.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_svoc    =  300.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_svoc    =  300.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_svoc    =  300.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_svoc    =  300.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_svoc   =  300.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_svoc   =  300.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_svoc   =  300.0_rk     ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_svoc   =  300.0_rk     ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_svoc   =  300.0_rk     ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_svoc   =  300.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Stress VOCs as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_svoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: agro_svoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: amat_svoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_svoc  = 1.0_rk

! Plant-Dependent emissions capacity/factors (EFs) for Other VOCs (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
! Other VOC (49 compounds): See Table 1 of Guenther et al. (2012)
        REAL(rk),          PARAMETER     :: ef1_ovoc    =  140.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef2_ovoc    =  140.0_rk     ! Needleleaf Evergreen Boreal Tree
        REAL(rk),          PARAMETER     :: ef3_ovoc    =  140.0_rk     ! Needleleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef4_ovoc    =  140.0_rk     ! Broadleaf Evergreen Tropical Tree
        REAL(rk),          PARAMETER     :: ef5_ovoc    =  140.0_rk     ! Broadleaf Evergreen Temperate Tree
        REAL(rk),          PARAMETER     :: ef6_ovoc    =  140.0_rk     ! Broadleaf Deciduous Tropical Tree
        REAL(rk),          PARAMETER     :: ef7_ovoc    =  140.0_rk     ! Broadleaf Deciduous Temperate Tree
        REAL(rk),          PARAMETER     :: ef8_ovoc    =  140.0_rk     ! Broadleaf Deciduous Boreal Tree
        REAL(rk),          PARAMETER     :: ef9_ovoc    =  140.0_rk     ! Broadleaf Evergreen Temperate Shrub
        REAL(rk),          PARAMETER     :: ef10_ovoc   =  140.0_rk     ! Broadleaf Deciduous Temperate Shrub
        REAL(rk),          PARAMETER     :: ef11_ovoc   =  140.0_rk     ! Broadleaf Deciduous Boreal Shrub
        REAL(rk),          PARAMETER     :: ef12_ovoc   =  140.0_rk     ! Arctic C3 Grass
        REAL(rk),          PARAMETER     :: ef13_ovoc   =  140.0_rk     ! Cool C3 Grass
        REAL(rk),          PARAMETER     :: ef14_ovoc   =  140.0_rk     ! Warm C4 Grass
        REAL(rk),          PARAMETER     :: ef15_ovoc   =  140.0_rk     ! Crop1

!Empirical factors or coefficients for: growing, mature, and old/senescing foliage, for Other VOCs as per Table 4 of Guenther et al., 2012
        REAL(rk),          PARAMETER     :: anew_ovoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: agro_ovoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: amat_ovoc  = 1.0_rk
        REAL(rk),          PARAMETER     :: aold_ovoc  = 1.0_rk

! Species-Dependent Parameterized Canopy Model Parameters (Table 4 of Guenther et al., 2012)
        REAL(rk),          PARAMETER     :: ldf_isop         =  1.0_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_isop        =  0.13_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_isop         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_isop         =  2.0_rk     !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_isop         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_isop         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_isop        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_isop         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_isop         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_isop        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_isop         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_isop         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_isop        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_isop         =  1.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_isop         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_isop        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_myrc         =  0.6_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_myrc        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_myrc         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_myrc         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_myrc         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_myrc         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_myrc        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_myrc         =  5.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_myrc         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_myrc        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_myrc         =  5.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_myrc         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_myrc        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_myrc         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_myrc         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_myrc        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_sabi         =  0.6_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_sabi        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_sabi         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_sabi         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_sabi         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_sabi         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_sabi        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_sabi         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_sabi         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_sabi        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_sabi         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_sabi         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_sabi        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_sabi         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_sabi         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_sabi        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_limo         =  0.2_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_limo        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_limo         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_limo         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_limo         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_limo         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_limo        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_limo         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_limo         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_limo        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_limo         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_limo         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_limo        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_limo         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_limo         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_limo        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_care         =  0.2_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_care        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_care         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_care         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_care         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_care         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_care        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_care         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_care         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_care        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_care         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_care         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_care        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_care         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_care         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_care        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_ocim         =  0.8_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_ocim        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_ocim         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_ocim         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_ocim         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_ocim         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_ocim        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_ocim         =  5.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_ocim         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_ocim        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_ocim         =  5.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_ocim         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_ocim        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_ocim         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_ocim         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_ocim        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_bpin         =  0.2_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_bpin        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_bpin         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_bpin         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_bpin         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_bpin         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_bpin        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_bpin         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_bpin         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_bpin        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_bpin         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_bpin         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_bpin        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_bpin         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_bpin         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_bpin        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_apin         =  0.6_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_apin        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_apin         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_apin         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_apin         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_apin         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_apin        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_apin         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_apin         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_apin        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_apin         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_apin         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_apin        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_apin         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_apin         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_apin        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_mono         =  0.4_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_mono        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_mono         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_mono         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_mono         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_mono         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_mono        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_mono         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_mono         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_mono        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_mono         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_mono         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_mono        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_mono         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_mono         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_mono        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_farn         =  0.5_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_farn        =  0.17_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_farn         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_farn         =  2.37_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_farn         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_farn         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_farn        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_farn         =  5.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_farn         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_farn        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_farn         =  5.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_farn         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_farn        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_farn         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_farn         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_farn        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_cary         =  0.5_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_cary        =  0.17_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_cary         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_cary         =  2.37_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_cary         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_cary         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_cary        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_cary         =  5.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_cary         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_cary        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_cary         =  5.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_cary         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_cary        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_cary         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_cary         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_cary        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_sesq         =  0.5_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_sesq        =  0.17_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_sesq         =  130.0_rk   !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_sesq         =  2.37_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_sesq         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_sesq         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_sesq        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_sesq         =  5.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_sesq         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_sesq        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_sesq         =  5.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_sesq         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_sesq        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_sesq         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_sesq         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_sesq        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_mbol         =  1.0_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_mbol        =  0.13_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_mbol         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_mbol         =  2.0_rk     !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_mbol         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_mbol         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_mbol        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_mbol         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_mbol         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_mbol        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_mbol         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_mbol         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_mbol        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_mbol         =  1.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_mbol         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_mbol        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_meth         =  0.8_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_meth        =  0.08_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_meth         =  60.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_meth         =  1.6_rk     !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_meth         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_meth         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_meth        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_meth         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_meth         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_meth        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_meth         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_meth         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_meth        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_meth         =  1.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_meth         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_meth        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_acet         =  0.2_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_acet        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_acet         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_acet         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_acet         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_acet         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_acet        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_acet         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_acet         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_acet        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_acet         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_acet         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_acet        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_acet         =  1.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_acet         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_acet        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_co           =  1.0_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_co          =  0.08_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_co           =  60.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_co           =  1.6_rk     !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_co           =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_co           =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_co          =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_co           =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_co           =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_co          =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_co           =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_co           =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_co          =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_co           =  1.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_co           =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_co          =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_bvoc         =  0.8_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_bvoc        =  0.13_rk    !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_bvoc         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_bvoc         =  2.0_rk     !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_bvoc         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_bvoc         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_bvoc        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_bvoc         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_bvoc         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_bvoc        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_bvoc         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_bvoc         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_bvoc        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_bvoc         =  1.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_bvoc         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_bvoc        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_svoc         =  0.8_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_svoc        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_svoc         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_svoc         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_svoc         =  5.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_svoc         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_svoc        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_svoc         =  5.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_svoc         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_svoc        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_svoc         =  5.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_svoc         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_svoc        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_svoc         =  5.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_svoc         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_svoc        =  8.0_rk     !delta threshold for high wind stress (m/s)

        REAL(rk),          PARAMETER     :: ldf_ovoc         =  0.2_rk     !Light-dependent fraction
        REAL(rk),          PARAMETER     :: beta_ovoc        =  0.1_rk     !Empirical coefficient for temperature dependence of light-independent fraction
        REAL(rk),          PARAMETER     :: ct1_ovoc         =  80.0_rk    !Activation energy (kJ/mol)
        REAL(rk),          PARAMETER     :: ceo_ovoc         =  1.83_rk    !Empirical coefficient
        REAL(rk),          PARAMETER     :: caq_ovoc         =  1.0_rk     !coefficient for poor Air Quality stress
        REAL(rk),          PARAMETER     :: taq_ovoc         =  20.0_rk    !threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: dtaq_ovoc        =  30.0_rk    !delta threshold for poor Air Quality stress (ppm-hours)
        REAL(rk),          PARAMETER     :: cht_ovoc         =  1.0_rk     !coefficient for high temperature stress
        REAL(rk),          PARAMETER     :: tht_ovoc         =  313.15_rk  !threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: dtht_ovoc        =  8.0_rk     !delta threshold for high temperature stress (K)
        REAL(rk),          PARAMETER     :: clt_ovoc         =  1.0_rk     !coefficient for low temperature stress
        REAL(rk),          PARAMETER     :: tlt_ovoc         =  283.15_rk  !threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: dtlt_ovoc        =  8.0_rk     !delta threshold for low temperature stress (K)
        REAL(rk),          PARAMETER     :: chw_ovoc         =  1.0_rk     !coefficient for high wind stress
        REAL(rk),          PARAMETER     :: thw_ovoc         =  12.0_rk    !threshold for high wind stress (m/s)
        REAL(rk),          PARAMETER     :: dthw_ovoc        =  8.0_rk     !delta threshold for high wind stress (m/s)

! Set tree and species dependent coefficients
        if (emi_ind .eq. 1 ) then
            ldf  = ldf_isop
            beta = beta_isop
            ct1  = ct1_isop
            ceo  = ceo_isop
            ef1  = ef1_isop
            ef2  = ef2_isop
            ef3  = ef3_isop
            ef4  = ef4_isop
            ef5  = ef5_isop
            ef6  = ef6_isop
            ef7  = ef7_isop
            ef8  = ef8_isop
            ef9  = ef9_isop
            ef10 = ef10_isop
            ef11 = ef11_isop
            ef12 = ef12_isop
            ef13 = ef13_isop
            ef14 = ef14_isop
            ef15 = ef15_isop
            anew  = anew_isop
            agro  = agro_isop
            amat  = amat_isop
            aold  = aold_isop
            caq   = caq_isop
            taq   = taq_isop
            dtaq  = dtaq_isop
            cht   = cht_isop
            tht   = tht_isop
            dtht  = dtht_isop
            clt   = clt_isop
            tlt   = tlt_isop
            dtlt  = dtlt_isop
            chw   = chw_isop
            thw   = thw_isop
            dthw  = dthw_isop
        else if (emi_ind .eq. 2 ) then
            ldf  = ldf_myrc
            beta = beta_myrc
            ct1 = ct1_myrc
            ceo = ceo_myrc
            ef1  = ef1_myrc
            ef2  = ef2_myrc
            ef3  = ef3_myrc
            ef4  = ef4_myrc
            ef5  = ef5_myrc
            ef6  = ef6_myrc
            ef7  = ef7_myrc
            ef8  = ef8_myrc
            ef9  = ef9_myrc
            ef10 = ef10_myrc
            ef11 = ef11_myrc
            ef12 = ef12_myrc
            ef13 = ef13_myrc
            ef14 = ef14_myrc
            ef15 = ef15_myrc
            anew  = anew_myrc
            agro  = agro_myrc
            amat  = amat_myrc
            aold  = aold_myrc
            caq   = caq_myrc
            taq   = taq_myrc
            dtaq  = dtaq_myrc
            cht   = cht_myrc
            tht   = tht_myrc
            dtht  = dtht_myrc
            clt   = clt_myrc
            tlt   = tlt_myrc
            dtlt  = dtlt_myrc
            chw   = chw_myrc
            thw   = thw_myrc
            dthw  = dthw_myrc
        else if (emi_ind .eq. 3 ) then
            ldf  = ldf_sabi
            beta = beta_sabi
            ct1 = ct1_sabi
            ceo = ceo_sabi
            ef1  = ef1_sabi
            ef2  = ef2_sabi
            ef3  = ef3_sabi
            ef4  = ef4_sabi
            ef5  = ef5_sabi
            ef6  = ef6_sabi
            ef7  = ef7_sabi
            ef8  = ef8_sabi
            ef9  = ef9_sabi
            ef10 = ef10_sabi
            ef11 = ef11_sabi
            ef12 = ef12_sabi
            ef13 = ef13_sabi
            ef14 = ef14_sabi
            ef15 = ef15_sabi
            anew  = anew_sabi
            agro  = agro_sabi
            amat  = amat_sabi
            aold  = aold_sabi
            caq   = caq_sabi
            taq   = taq_sabi
            dtaq  = dtaq_sabi
            cht   = cht_sabi
            tht   = tht_sabi
            dtht  = dtht_sabi
            clt   = clt_sabi
            tlt   = tlt_sabi
            dtlt  = dtlt_sabi
            chw   = chw_sabi
            thw   = thw_sabi
            dthw  = dthw_sabi
        else if (emi_ind .eq. 4 ) then
            ldf  = ldf_limo
            beta = beta_limo
            ct1 = ct1_limo
            ceo = ceo_limo
            ef1  = ef1_limo
            ef2  = ef2_limo
            ef3  = ef3_limo
            ef4  = ef4_limo
            ef5  = ef5_limo
            ef6  = ef6_limo
            ef7  = ef7_limo
            ef8  = ef8_limo
            ef9  = ef9_limo
            ef10 = ef10_limo
            ef11 = ef11_limo
            ef12 = ef12_limo
            ef13 = ef13_limo
            ef14 = ef14_limo
            ef15 = ef15_limo
            anew  = anew_limo
            agro  = agro_limo
            amat  = amat_limo
            aold  = aold_limo
            caq   = caq_limo
            taq   = taq_limo
            dtaq  = dtaq_limo
            cht   = cht_limo
            tht   = tht_limo
            dtht  = dtht_limo
            clt   = clt_limo
            tlt   = tlt_limo
            dtlt  = dtlt_limo
            chw   = chw_limo
            thw   = thw_limo
            dthw  = dthw_limo
        else if (emi_ind .eq. 5 ) then
            ldf  = ldf_care
            beta = beta_care
            ct1 = ct1_care
            ceo = ceo_care
            ef1  = ef1_care
            ef2  = ef2_care
            ef3  = ef3_care
            ef4  = ef4_care
            ef5  = ef5_care
            ef6  = ef6_care
            ef7  = ef7_care
            ef8  = ef8_care
            ef9  = ef9_care
            ef10 = ef10_care
            ef11 = ef11_care
            ef12 = ef12_care
            ef13 = ef13_care
            ef14 = ef14_care
            ef15 = ef15_care
            anew  = anew_care
            agro  = agro_care
            amat  = amat_care
            aold  = aold_care
            caq   = caq_care
            taq   = taq_care
            dtaq  = dtaq_care
            cht   = cht_care
            tht   = tht_care
            dtht  = dtht_care
            clt   = clt_care
            tlt   = tlt_care
            dtlt  = dtlt_care
            chw   = chw_care
            thw   = thw_care
            dthw  = dthw_care
        else if (emi_ind .eq. 6 ) then
            ldf  = ldf_ocim
            beta = beta_ocim
            ct1 = ct1_ocim
            ceo = ceo_ocim
            ef1  = ef1_ocim
            ef2  = ef2_ocim
            ef3  = ef3_ocim
            ef4  = ef4_ocim
            ef5  = ef5_ocim
            ef6  = ef6_ocim
            ef7  = ef7_ocim
            ef8  = ef8_ocim
            ef9  = ef9_ocim
            ef10 = ef10_ocim
            ef11 = ef11_ocim
            ef12 = ef12_ocim
            ef13 = ef13_ocim
            ef14 = ef14_ocim
            ef15 = ef15_ocim
            anew  = anew_ocim
            agro  = agro_ocim
            amat  = amat_ocim
            aold  = aold_ocim
            caq   = caq_ocim
            taq   = taq_ocim
            dtaq  = dtaq_ocim
            cht   = cht_ocim
            tht   = tht_ocim
            dtht  = dtht_ocim
            clt   = clt_ocim
            tlt   = tlt_ocim
            dtlt  = dtlt_ocim
            chw   = chw_ocim
            thw   = thw_ocim
            dthw  = dthw_ocim
        else if (emi_ind .eq. 7 ) then
            ldf  = ldf_bpin
            beta = beta_bpin
            ct1 = ct1_bpin
            ceo = ceo_bpin
            ef1  = ef1_bpin
            ef2  = ef2_bpin
            ef3  = ef3_bpin
            ef4  = ef4_bpin
            ef5  = ef5_bpin
            ef6  = ef6_bpin
            ef7  = ef7_bpin
            ef8  = ef8_bpin
            ef9  = ef9_bpin
            ef10 = ef10_bpin
            ef11 = ef11_bpin
            ef12 = ef12_bpin
            ef13 = ef13_bpin
            ef14 = ef14_bpin
            ef15 = ef15_bpin
            anew  = anew_bpin
            agro  = agro_bpin
            amat  = amat_bpin
            aold  = aold_bpin
            caq   = caq_bpin
            taq   = taq_bpin
            dtaq  = dtaq_bpin
            cht   = cht_bpin
            tht   = tht_bpin
            dtht  = dtht_bpin
            clt   = clt_bpin
            tlt   = tlt_bpin
            dtlt  = dtlt_bpin
            chw   = chw_bpin
            thw   = thw_bpin
            dthw  = dthw_bpin
        else if (emi_ind .eq. 8 ) then
            ldf  = ldf_apin
            beta = beta_apin
            ct1 = ct1_apin
            ceo = ceo_apin
            ef1  = ef1_apin
            ef2  = ef2_apin
            ef3  = ef3_apin
            ef4  = ef4_apin
            ef5  = ef5_apin
            ef6  = ef6_apin
            ef7  = ef7_apin
            ef8  = ef8_apin
            ef9  = ef9_apin
            ef10 = ef10_apin
            ef11 = ef11_apin
            ef12 = ef12_apin
            ef13 = ef13_apin
            ef14 = ef14_apin
            ef15 = ef15_apin
            anew  = anew_apin
            agro  = agro_apin
            amat  = amat_apin
            aold  = aold_apin
            caq   = caq_apin
            taq   = taq_apin
            dtaq  = dtaq_apin
            cht   = cht_apin
            tht   = tht_apin
            dtht  = dtht_apin
            clt   = clt_apin
            tlt   = tlt_apin
            dtlt  = dtlt_apin
            chw   = chw_apin
            thw   = thw_apin
            dthw  = dthw_apin
        else if (emi_ind .eq. 9 ) then
            ldf  = ldf_mono
            beta = beta_mono
            ct1 = ct1_mono
            ceo = ceo_mono
            ef1  = ef1_mono
            ef2  = ef2_mono
            ef3  = ef3_mono
            ef4  = ef4_mono
            ef5  = ef5_mono
            ef6  = ef6_mono
            ef7  = ef7_mono
            ef8  = ef8_mono
            ef9  = ef9_mono
            ef10 = ef10_mono
            ef11 = ef11_mono
            ef12 = ef12_mono
            ef13 = ef13_mono
            ef14 = ef14_mono
            ef15 = ef15_mono
            anew  = anew_mono
            agro  = agro_mono
            amat  = amat_mono
            aold  = aold_mono
            caq   = caq_mono
            taq   = taq_mono
            dtaq  = dtaq_mono
            cht   = cht_mono
            tht   = tht_mono
            dtht  = dtht_mono
            clt   = clt_mono
            tlt   = tlt_mono
            dtlt  = dtlt_mono
            chw   = chw_mono
            thw   = thw_mono
            dthw  = dthw_mono
        else if (emi_ind .eq. 10 ) then
            ldf  = ldf_farn
            beta = beta_farn
            ct1 = ct1_farn
            ceo = ceo_farn
            ef1  = ef1_farn
            ef2  = ef2_farn
            ef3  = ef3_farn
            ef4  = ef4_farn
            ef5  = ef5_farn
            ef6  = ef6_farn
            ef7  = ef7_farn
            ef8  = ef8_farn
            ef9  = ef9_farn
            ef10 = ef10_farn
            ef11 = ef11_farn
            ef12 = ef12_farn
            ef13 = ef13_farn
            ef14 = ef14_farn
            ef15 = ef15_farn
            anew  = anew_farn
            agro  = agro_farn
            amat  = amat_farn
            aold  = aold_farn
            caq   = caq_farn
            taq   = taq_farn
            dtaq  = dtaq_farn
            cht   = cht_farn
            tht   = tht_farn
            dtht  = dtht_farn
            clt   = clt_farn
            tlt   = tlt_farn
            dtlt  = dtlt_farn
            chw   = chw_farn
            thw   = thw_farn
            dthw  = dthw_farn
        else if (emi_ind .eq. 11 ) then
            ldf  = ldf_cary
            beta = beta_cary
            ct1 = ct1_cary
            ceo = ceo_cary
            ef1  = ef1_cary
            ef2  = ef2_cary
            ef3  = ef3_cary
            ef4  = ef4_cary
            ef5  = ef5_cary
            ef6  = ef6_cary
            ef7  = ef7_cary
            ef8  = ef8_cary
            ef9  = ef9_cary
            ef10 = ef10_cary
            ef11 = ef11_cary
            ef12 = ef12_cary
            ef13 = ef13_cary
            ef14 = ef14_cary
            ef15 = ef15_cary
            anew  = anew_cary
            agro  = agro_cary
            amat  = amat_cary
            aold  = aold_cary
            caq   = caq_cary
            taq   = taq_cary
            dtaq  = dtaq_cary
            cht   = cht_cary
            tht   = tht_cary
            dtht  = dtht_cary
            clt   = clt_cary
            tlt   = tlt_cary
            dtlt  = dtlt_cary
            chw   = chw_cary
            thw   = thw_cary
            dthw  = dthw_cary
        else if (emi_ind .eq. 12 ) then
            ldf  = ldf_sesq
            beta = beta_sesq
            ct1 = ct1_sesq
            ceo = ceo_sesq
            ef1  = ef1_sesq
            ef2  = ef2_sesq
            ef3  = ef3_sesq
            ef4  = ef4_sesq
            ef5  = ef5_sesq
            ef6  = ef6_sesq
            ef7  = ef7_sesq
            ef8  = ef8_sesq
            ef9  = ef9_sesq
            ef10 = ef10_sesq
            ef11 = ef11_sesq
            ef12 = ef12_sesq
            ef13 = ef13_sesq
            ef14 = ef14_sesq
            ef15 = ef15_sesq
            anew  = anew_sesq
            agro  = agro_sesq
            amat  = amat_sesq
            aold  = aold_sesq
            caq   = caq_sesq
            taq   = taq_sesq
            dtaq  = dtaq_sesq
            cht   = cht_sesq
            tht   = tht_sesq
            dtht  = dtht_sesq
            clt   = clt_sesq
            tlt   = tlt_sesq
            dtlt  = dtlt_sesq
            chw   = chw_sesq
            thw   = thw_sesq
            dthw  = dthw_sesq
        else if (emi_ind .eq. 13 ) then
            ldf  = ldf_mbol
            beta = beta_mbol
            ct1 = ct1_mbol
            ceo = ceo_mbol
            ef1  = ef1_mbol
            ef2  = ef2_mbol
            ef3  = ef3_mbol
            ef4  = ef4_mbol
            ef5  = ef5_mbol
            ef6  = ef6_mbol
            ef7  = ef7_mbol
            ef8  = ef8_mbol
            ef9  = ef9_mbol
            ef10 = ef10_mbol
            ef11 = ef11_mbol
            ef12 = ef12_mbol
            ef13 = ef13_mbol
            ef14 = ef14_mbol
            ef15 = ef15_mbol
            anew  = anew_mbol
            agro  = agro_mbol
            amat  = amat_mbol
            aold  = aold_mbol
            caq   = caq_mbol
            taq   = taq_mbol
            dtaq  = dtaq_mbol
            cht   = cht_mbol
            tht   = tht_mbol
            dtht  = dtht_mbol
            clt   = clt_mbol
            tlt   = tlt_mbol
            dtlt  = dtlt_mbol
            chw   = chw_mbol
            thw   = thw_mbol
            dthw  = dthw_mbol
        else if (emi_ind .eq. 14 ) then
            ldf  = ldf_meth
            beta = beta_meth
            ct1 = ct1_meth
            ceo = ceo_meth
            ef1  = ef1_meth
            ef2  = ef2_meth
            ef3  = ef3_meth
            ef4  = ef4_meth
            ef5  = ef5_meth
            ef6  = ef6_meth
            ef7  = ef7_meth
            ef8  = ef8_meth
            ef9  = ef9_meth
            ef10 = ef10_meth
            ef11 = ef11_meth
            ef12 = ef12_meth
            ef13 = ef13_meth
            ef14 = ef14_meth
            ef15 = ef15_meth
            anew  = anew_meth
            agro  = agro_meth
            amat  = amat_meth
            aold  = aold_meth
            caq   = caq_meth
            taq   = taq_meth
            dtaq  = dtaq_meth
            cht   = cht_meth
            tht   = tht_meth
            dtht  = dtht_meth
            clt   = clt_meth
            tlt   = tlt_meth
            dtlt  = dtlt_meth
            chw   = chw_meth
            thw   = thw_meth
            dthw  = dthw_meth
        else if (emi_ind .eq. 15 ) then
            ldf  = ldf_acet
            beta = beta_acet
            ct1 = ct1_acet
            ceo = ceo_acet
            ef1  = ef1_acet
            ef2  = ef2_acet
            ef3  = ef3_acet
            ef4  = ef4_acet
            ef5  = ef5_acet
            ef6  = ef6_acet
            ef7  = ef7_acet
            ef8  = ef8_acet
            ef9  = ef9_acet
            ef10 = ef10_acet
            ef11 = ef11_acet
            ef12 = ef12_acet
            ef13 = ef13_acet
            ef14 = ef14_acet
            ef15 = ef15_acet
            anew  = anew_acet
            agro  = agro_acet
            amat  = amat_acet
            aold  = aold_acet
            caq   = caq_acet
            taq   = taq_acet
            dtaq  = dtaq_acet
            cht   = cht_acet
            tht   = tht_acet
            dtht  = dtht_acet
            clt   = clt_acet
            tlt   = tlt_acet
            dtlt  = dtlt_acet
            chw   = chw_acet
            thw   = thw_acet
            dthw  = dthw_acet
        else if (emi_ind .eq. 16 ) then
            ldf  = ldf_co
            beta = beta_co
            ct1 = ct1_co
            ceo = ceo_co
            ef1  = ef1_co
            ef2  = ef2_co
            ef3  = ef3_co
            ef4  = ef4_co
            ef5  = ef5_co
            ef6  = ef6_co
            ef7  = ef7_co
            ef8  = ef8_co
            ef9  = ef9_co
            ef10 = ef10_co
            ef11 = ef11_co
            ef12 = ef12_co
            ef13 = ef13_co
            ef14 = ef14_co
            ef15 = ef15_co
            anew  = anew_co
            agro  = agro_co
            amat  = amat_co
            aold  = aold_co
            caq   = caq_co
            taq   = taq_co
            dtaq  = dtaq_co
            cht   = cht_co
            tht   = tht_co
            dtht  = dtht_co
            clt   = clt_co
            tlt   = tlt_co
            dtlt  = dtlt_co
            chw   = chw_co
            thw   = thw_co
            dthw  = dthw_co
        else if (emi_ind .eq. 17 ) then
            ldf  = ldf_bvoc
            beta = beta_bvoc
            ct1 = ct1_bvoc
            ceo = ceo_bvoc
            ef1  = ef1_bvoc
            ef2  = ef2_bvoc
            ef3  = ef3_bvoc
            ef4  = ef4_bvoc
            ef5  = ef5_bvoc
            ef6  = ef6_bvoc
            ef7  = ef7_bvoc
            ef8  = ef8_bvoc
            ef9  = ef9_bvoc
            ef10 = ef10_bvoc
            ef11 = ef11_bvoc
            ef12 = ef12_bvoc
            ef13 = ef13_bvoc
            ef14 = ef14_bvoc
            ef15 = ef15_bvoc
            anew  = anew_bvoc
            agro  = agro_bvoc
            amat  = amat_bvoc
            aold  = aold_bvoc
            caq   = caq_bvoc
            taq   = taq_bvoc
            dtaq  = dtaq_bvoc
            cht   = cht_bvoc
            tht   = tht_bvoc
            dtht  = dtht_bvoc
            clt   = clt_bvoc
            tlt   = tlt_bvoc
            dtlt  = dtlt_bvoc
            chw   = chw_bvoc
            thw   = thw_bvoc
            dthw  = dthw_bvoc
        else if (emi_ind .eq. 18 ) then
            ldf  = ldf_svoc
            beta = beta_svoc
            ct1 = ct1_svoc
            ceo = ceo_svoc
            ef1  = ef1_svoc
            ef2  = ef2_svoc
            ef3  = ef3_svoc
            ef4  = ef4_svoc
            ef5  = ef5_svoc
            ef6  = ef6_svoc
            ef7  = ef7_svoc
            ef8  = ef8_svoc
            ef9  = ef9_svoc
            ef10 = ef10_svoc
            ef11 = ef11_svoc
            ef12 = ef12_svoc
            ef13 = ef13_svoc
            ef14 = ef14_svoc
            ef15 = ef15_svoc
            anew  = anew_svoc
            agro  = agro_svoc
            amat  = amat_svoc
            aold  = aold_svoc
            caq   = caq_svoc
            taq   = taq_svoc
            dtaq  = dtaq_svoc
            cht   = cht_svoc
            tht   = tht_svoc
            dtht  = dtht_svoc
            clt   = clt_svoc
            tlt   = tlt_svoc
            dtlt  = dtlt_svoc
            chw   = chw_svoc
            thw   = thw_svoc
            dthw  = dthw_svoc
        else   ! EMI_IND = 19
            ldf  = ldf_ovoc
            beta = beta_ovoc
            ct1 = ct1_ovoc
            ceo = ceo_ovoc
            ef1  = ef1_ovoc
            ef2  = ef2_ovoc
            ef3  = ef3_ovoc
            ef4  = ef4_ovoc
            ef5  = ef5_ovoc
            ef6  = ef6_ovoc
            ef7  = ef7_ovoc
            ef8  = ef8_ovoc
            ef9  = ef9_ovoc
            ef10 = ef10_ovoc
            ef11 = ef11_ovoc
            ef12 = ef12_ovoc
            ef13 = ef13_ovoc
            ef14 = ef14_ovoc
            ef15 = ef15_ovoc
            anew  = anew_ovoc
            agro  = agro_ovoc
            amat  = amat_ovoc
            aold  = aold_ovoc
            caq   = caq_ovoc
            taq   = taq_ovoc
            dtaq  = dtaq_ovoc
            cht   = cht_ovoc
            tht   = tht_ovoc
            dtht  = dtht_ovoc
            clt   = clt_ovoc
            tlt   = tlt_ovoc
            dtlt  = dtlt_ovoc
            chw   = chw_ovoc
            thw   = thw_ovoc
            dthw  = dthw_ovoc
        end if

        if (lu_opt .eq. 0 .or. lu_opt .eq. 1) then !VIIRS or MODIS  LU types

! Simple MEGAN (Table 3 in Guenther et al., 2012) PFT to VIIRS/MODIS VTYPE mapping
            if (vtype .eq. 1) then !VIIRS Cat 1 Evergreen Needleleaf
                !--> Average Needleleaf Evergreen Temperate Tree and Needleleaf Evergreen Boreal Tree

                ef = (ef1+ef2)/2.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 6.706_rk
                rootb = 2.175_rk

            else if (vtype .eq. 2) then !VIIRS/MODIS Cat 2 Evergreen Broadleaf
                !--> Average Broadleaf Evergreen Tropical Tree and Broadleaf Evergreen Temperate Tree

                ef = (ef4+ef5)/2.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 7.344_rk
                rootb = 1.303_rk

            else if (vtype .eq. 3) then !VIIRS/MODIS Cat 3 Deciduous Needleaf
                !--> Average Needleleaf Deciduous Boreal Tree

                ef = ef3
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 7.066_rk
                rootb = 1.953_rk

            else if (vtype .eq. 4) then !VIIRS/MODIS Cat 4 Deciduous Broadleaf
                !--> Average Broadleaf Deciduous Tropical Tree, Broadleaf Deciduous Temperate Tree,
                ! and Broadleaf Deciduous Boreal Tree

                ef = (ef6+ef7+ef8)/3.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 5.990_rk
                rootb = 1.955_rk

            else if (vtype .eq. 5) then !VIIRS/MODIS Cat 5 Mixed forests
                !--> Avearge of all above EF1-EF8 PFTs.

                ef = (ef1+ef2+ef3+ef4+ef5+ef6+ef7+ef8)/8.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 4.453_rk
                rootb = 1.631_rk

            else if (vtype .ge. 6 .and. vtype .le. 7) then !VIIRS/MODIS Cat 6-7 Closed/Open Shrublands
                !--> Avearge Broadleaf Evergreen Temperate Shrub, Broadleaf Deciduous Temperate Shrub,
                ! and Broadleaf Deciduous Boreal Shrub

                ef = (ef9+ef10+ef11)/3.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = (6.326_rk + 7.718_rk)/2.0_rk
                rootb = (1.567_rk + 1.262_rk)/2.0_rk

            else if (vtype .ge. 8 .and. vtype .le. 11) then !VIIRS/MODIS Cat 8-10 Savannas and Grasslands
                !--> Avearge Arctic C3 Grass, Cool C3 Grass, Warm C4 Grass)

                ef = (ef12+ef13+ef14)/3.0_rk
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = (7.604_rk + 8.235_rk + 10.74_rk)/3.0_rk
                rootb = (2.300_rk + 1.627_rk + 2.608_rk)/3.0_rk

            else if (vtype .eq. 12 ) then !VIIRS/MODIS Cat 12 Croplands
                !--> Crop1

                ef = ef15
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 5.558_rk
                rootb = 2.614_rk

            else if (vtype .eq. 14 ) then !VIIRS/MODIS Cat 14 Croplands/Natural Mosaic
                !--> Crop1

                ef = ef15
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = 5.558_rk
                rootb = 2.614_rk

            else if (vtype .ge. 18 .and. vtype .le. 19) then !VIIRS/MODIS Cat 18-19 Wooded/Mixed Tundra
                !--> Crop1

                ef = (ef9+ef10+ef11)/3.0_rk  !Assume EF for wooded/mixed tundra = shrublands above
                !Set PFT dependent A and B coefficients for cumulative root depth fraction (Zeng 2001)
                !See Table 2 for IGPB Classification at:  https://journals.ametsoc.org/view/journals/hydr/2/5/1525-7541_2001_002_0525_gvrdfl_2_0_co_2.xml
                roota = (6.326_rk + 7.718_rk)/2.0_rk
                rootb = (1.567_rk + 1.262_rk)/2.0_rk

            else

                ef = 0.0_rk  !set EF=0 for non-valid VTYPEs

            end if

        else
            write(*,*)  'Wrong LU_OPT choice of ', lu_opt, 'in namelist, only VIIRS (0) or MODIS (1) available right now...exiting'
            call exit(2)
        end if

    END SUBROUTINE canopy_biop


end module canopy_bioparm_mod
```


