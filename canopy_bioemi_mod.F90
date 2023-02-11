module canopy_bioemi_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_BIO( ZK, FCLAI, FCH, LAI, CLU, COSZEN, SFCRAD, &
                          TEMP2, LU_OPT, VTYPE, EMI_OUT)
!                            , EMI_NAME, EMI_OUT ) !TBD

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
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk,rgasuniv   !constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal
        use canopy_phot_mod

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: ZK(:)           ! Input model heights (m)
        REAL(RK),    INTENT( IN )  :: FCLAI(:)        ! Input Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )  :: FCH             ! Model input canopy height (m)
        REAL(RK),    INTENT( IN )  :: LAI             ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )  :: CLU             ! Model input Clumping Index
        REAL(RK),    INTENT( IN )  :: COSZEN          ! Model input Cosine Solar Zenith Angle
        REAL(RK),    INTENT( IN )  :: SFCRAD          ! Model input Instantaneous surface downward shortwave flux (W/m2)
        REAL(RK),    INTENT( IN )  :: TEMP2           ! Model input 2-m Temperature (K)
        INTEGER,     INTENT( IN )  :: LU_OPT          ! integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
        INTEGER,     INTENT( IN )  :: VTYPE           ! Grid cell dominant vegetation type
        REAL(RK),    INTENT( OUT ) :: EMI_OUT(:)      ! Output emissions (kg /m2 s)
!        REAL(RK)    :: EMI_OUT(SIZE(ZK)) !test bioemis

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
        REAL(RK) :: EF                             ! Emission factor (EF) (ug/m2 hr)
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

! Plant-Dependent emissions capacity/factors (EFs) (Tables 2-3 of Guenther et al., 2012) (ug/m2 hr)
        REAL(RK),          PARAMETER     :: EF1    =  600.0_rk     ! Needleleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF2    =  3000.0_rk    ! Needleleaf Evergreen Boreal Tree
        REAL(RK),          PARAMETER     :: EF3    =  1.0_rk       ! Needleleaf Deciduous Boreal Tree
        REAL(RK),          PARAMETER     :: EF4    =  7000.0_rk    ! Broadleaf Evergreen Tropical Tree
        REAL(RK),          PARAMETER     :: EF5    =  10000.0_rk   ! Broadleaf Evergreen Temperate Tree
        REAL(RK),          PARAMETER     :: EF6    =  7000.0_rk    ! Broadleaf Deciduous Tropical Tree
        REAL(RK),          PARAMETER     :: EF7    =  10000.0_rk   ! Broadleaf Deciduous Temperate Tree
        REAL(RK),          PARAMETER     :: EF8    =  11000.0_rk   ! Broadleaf Deciduous Boreal Tree                                                                          
        REAL(RK),          PARAMETER     :: EF9    =  2000.0_rk    ! Broadleaf Evergreen Temperate Shrub
        REAL(RK),          PARAMETER     :: EF10   =  4000.0_rk    ! Broadleaf Deciduous Temperate Shrub
        REAL(RK),          PARAMETER     :: EF11   =  4000.0_rk    ! Broadleaf Deciduous Boreal Shrub
        REAL(RK),          PARAMETER     :: EF12   =  1600.0_rk    ! Arctic C3 Grass
        REAL(RK),          PARAMETER     :: EF13   =  800.0_rk     ! Cool C3 Grass
        REAL(RK),          PARAMETER     :: EF14   =  200.0_rk     ! Warm C4 Grass
        REAL(RK),          PARAMETER     :: EF15   =  1.0_rk       ! Crop1

! Species-Dependent Parameterized Canopy Model Parameters (Table 4 of Guenther et al., 2012)
        REAL(RK),          PARAMETER     :: CT1_ISOP         =  95.0_rk    !Activation energy (kJ/mol)
        REAL(RK),          PARAMETER     :: CEO_ISOP         =  2.0_rk     !Empirical coefficient
!        REAL(RK),          PARAMETER     :: CT1_MYRC         =  80.0_rk    !Activation energy (kJ/mol)
!        REAL(RK),          PARAMETER     :: CEO_MYRC         =  1.83_rk     !Empirical coefficient


!Calculate photolyis shading/correction factor through canopy, i.e., the fraction of sunlit leaves downward through canopy

       call canopy_phot(FCLAI, LAI, CLU, COSZEN, RJCF)

        ! Calculate PPFD sun and shade through canopy layers using photolysis reduction factors
        ! above canopy RJCF = 1 

       PPFD_SUN   = RJCF*SFCRAD*PPFD_CONST
       PPFD_SHADE = MAX((1.0-RJCF),0.0_rk)*SFCRAD*PPFD_CONST
!       print*, 'Sunlit = ', RJCF, 'Shaded=', MAX((1.0-RJCF),0.0_rk)
!       print*, 'PPFD_SHADE=',PPFD_SHADE
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

! Set species dependent coefficients

!TBD Check Emission Name
!       if (EMI_NAME .eq. "ISOP" ) then
          CT1 = CT1_ISOP
          CEO = CEO_ISOP
!       else if (EMI_NAME .eq. "MYRC" ) then 
!          CT1 = CT1_MYRC
!          CEO = CEO_MYRC
!       end if

       E_OPT = CEO * EXP(0.05_rk * (TLEAF24_AVE-297.0_rk)) * EXP(0.05_rk * (TLEAF240_AVE-297.0_rk))
!       print*, 'TLEAF_AVE = ', TLEAF_AVE, 'E_OPT = ', E_OPT, 'TLEAF_OPT = ', TLEAF_OPT

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
!       print*, 'PPFD240_SUN=',PPFD240_SUN,'PPFD240_SHADE=',PPFD240_SHADE
       GammaPPFD_SUN   = CP_SUN*((ALPHA_P_SUN*PPFD_SUN)/SQRT(1.0 + (ALPHA_P_SUN**2.0) * (PPFD_SUN**2.0)))
       GammaPPFD_SHADE = CP_SHADE*((ALPHA_P_SHADE*PPFD_SHADE)/SQRT(1.0 + (ALPHA_P_SHADE**2.0) * (PPFD_SHADE**2.0)))
       
       GammaPPFD_AVE = (GammaPPFD_SUN*RJCF) + (GammaPPFD_SHADE*(1.0-RJCF)) ! average = sum sun and shade weighted by sunlit fraction
!       print*, 'CP_SUN=',CP_SUN,'CP_SHADE=',CP_SHADE,'ALPHA_P_SUN=',ALPHA_P_SUN,'ALPHA_P_SHADE=',ALPHA_P_SHADE

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
           if (ZK(i) .gt. 0.0 .and. ZK(i) .le. FCH) then      ! above ground level and at/below canopy top
               FLAI(i) = (FCLAI(i+1) - FCLAI(i)) * LAI  !fractional LAI in layer (surrogate for using LAD to get emissions in layer)
               EMI_OUT(i) = FLAI(i) * EF * GammaTLEAF_AVE(i) * GammaPPFD_AVE(i)  ! (ug/m2 hr)
               EMI_OUT(i) = EMI_OUT(i) * 2.77778E-13 !TBD:  convert emissions output to kg/m2 s
           end if
!                print*, 'Z/Hc=',ZK(i)/FCH, 'Gamma*FLAI = ', GammaTLEAF_AVE(i)*GammaPPFD_AVE(i)* FLAI(i), &
!                        'VTYPE=', VTYPE, 'EF = ', EF, 'EMI_OUT (umol/m2 hr) =',EMI_OUT(i)/68.12_rk
!               print*, 'Z/Hc=',ZK(i)/FCH, 'FLAI = ', FLAI(i), 'EMI_OUT (ug/m2 hr) =',EMI_OUT(i)
       end do 
!               print*, SUM(EMI_OUT)/68.12_rk
!              print*, 'ZK/HCM=', ZK/FCH, 'GammaTLEAF_AVE=',GammaTLEAF_AVE, 'GammaPPFD_AVE=', GammaPPFD_AVE,'EMI_OUT =',EMI_OUT
!               print*,'EMI_OUT=',EMI_OUT

    END SUBROUTINE CANOPY_BIO

end module canopy_bioemi_mod
