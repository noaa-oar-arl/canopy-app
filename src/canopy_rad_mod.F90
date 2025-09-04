!> \file canopy_rad_mod.F90
!> \brief Canopy radiation and photosynthetic photon flux density calculations
!> \details This module contains routines for calculating sunlit/shaded fractions
!>          and photosynthetic photon flux density (PPFD) profiles within forest
!>          canopies using exponential decay models based on Silva et al. (2020).
!> \author P. C. Campbell
!> \date Jun 2023

!> \defgroup canopy_radiation Canopy Radiation Calculations
!> \brief Radiation attenuation and PPFD calculations within canopy

module canopy_rad_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !> \brief Calculate sunlit fraction using clumping index
    !> \details Computes sunlit/shaded fraction through canopy using photolysis
    !!          correction factor and clumping index
    !> \ingroup canopy_radiation
    !> \param[in] FCLAI Fractional cumulative LAI profile (dimensionless)
    !> \param[in] LAI Model input total Leaf Area Index (m²/m²)
    !> \param[in] CLU Model input Clumping Index (dimensionless)
    !> \param[in] COSZEN Model input Cosine Solar Zenith Angle (dimensionless)
    !> \param[out] FSUN Sunlit/Shaded fraction from photolysis correction factor
    !> \author P. C. Campbell
    !> \date Jun 2023
    !> \note Based on Bonan (2019) equation 14.18 for clumping correction
    SUBROUTINE CANOPY_FSUN_CLU( FCLAI, LAI, CLU, COSZEN, FSUN)

!-----------------------------------------------------------------------

! Description:
!     computes linear interpolation method for PPFD sun/shade in canopy.

! Preconditions:
!     in-canopy height, and model LAI, clumping index, and solar zenith angle

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/23 by PCC
!     Jun 2023 P.C. Campbell: Initial standalone PPFD linear subroutine based on
!                             Silva et al. (2020) exponential curve algorithms
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: RK     !constants for canopy models
        use canopy_phot_mod

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: FCLAI(:)          ! Input Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )       :: LAI               ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )       :: CLU               ! Model input Clumping Index
        REAL(RK),    INTENT( IN )       :: COSZEN            ! Model input Cosine Solar Zenith Angle
        REAL(RK),    INTENT( OUT )      :: FSUN(SIZE(FCLAI)) ! Sunlit/Shaded fraction from photolysis correction factor

!Calculate photolyis shading/correction factor through canopy, i.e., the fraction of sunlit leaves downward through canopy
!  `canopy_phot` gives relative direct beam irradiance,
!  which, multiplied by clumping index, gives sunlit fraction (e.g., Bonan 2019, eq. 14.18)

        call canopy_phot(FCLAI, LAI, CLU, COSZEN, FSUN)
        FSUN = FSUN * CLU

    END SUBROUTINE CANOPY_FSUN_CLU

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !> \brief Calculate PPFD profiles using exponential model
    !> \details Computes photosynthetic photon flux density for sunlit and shaded leaves
    !!          through canopy using exponential decay models from Silva et al. (2020)
    !> \ingroup canopy_radiation
    !> \param[in] ZK Input model heights (m)
    !> \param[in] FCH Model input canopy height (m)
    !> \param[in] SFCRAD Model input instantaneous surface downward shortwave flux (W/m²)
    !> \param[in] LAI Model input total Leaf Area Index (m²/m²)
    !> \param[in] FSUN Sunlit/Shaded fraction from photolysis correction factor
    !> \param[out] PPFD_SUN PPFD for sunlit leaves (μmol photons/m²/s)
    !> \param[out] PPFD_SHADE PPFD for shaded leaves (μmol photons/m²/s)
    !> \param[out] PPFD_AVE Average PPFD for sunlit and shaded leaves (μmol photons/m²/s)
    !> \author P. C. Campbell
    !> \date Jun 2023
    !> \note Based on Silva et al. (2020) 5-layer canopy exponential PPFD model
    !> \cite Silva, S. J., Heald, C. L., and Guenther, A. B. (2020). Development of a
    !!       reduced-complexity plant canopy physics surrogate model for use in chemical
    !!       transport models: a case study with GEOS-Chem v12.3.0. Geoscientific Model
    !!       Development, 13, 2569-2585. https://doi.org/10.5194/gmd-13-2569-2020
    SUBROUTINE CANOPY_PPFD_EXP( ZK, FCH, SFCRAD, LAI, FSUN, &
        PPFD_SUN, PPFD_SHADE, PPFD_AVE)

!-----------------------------------------------------------------------

! Description:
!     computes linear interpolation method for PPFD sun/shade in canopy.

! Preconditions:
!     in-canopy height, and model LAI, clumping index, and solar zenith angle

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/23 by PCC
!     Jun 2023 P.C. Campbell: Initial standalone PPFD linear subroutine based on
!                             Silva et al. (2020) exponential curve algorithms
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: RK     !constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: ZK(:)                          ! Input model heights (m)
        REAL(RK),    INTENT( IN )       :: FCH                            ! Model input canopy height (m)
        REAL(RK),    INTENT( IN )       :: SFCRAD                         ! Model input Instantaneous surface downward shortwave flux (W/m2)
        REAL(RK),    INTENT( IN )       :: LAI                            ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )       :: FSUN(:)                        ! Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( OUT )      :: PPFD_SUN(SIZE(ZK))             ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( OUT )      :: PPFD_SHADE(SIZE(ZK))           ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( OUT )      :: PPFD_AVE(SIZE(ZK))             ! Average PPFD for sunlit and shaded leaves (umol phot/m2 s)

!      LOCAL
        !> \brief Exponential PPFD regression coefficient C for sunlit leaves at level 1
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 1 (top of canopy)
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_1_SUN     =  1.083_rk

        !> \brief Exponential PPFD regression coefficient C for sunlit leaves at level 2
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 2
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_2_SUN     =  1.096_rk

        !> \brief Exponential PPFD regression coefficient C for sunlit leaves at level 3
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 3
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_3_SUN     =  1.104_rk

        !> \brief Exponential PPFD regression coefficient C for sunlit leaves at level 4
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 4
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_4_SUN     =  1.098_rk

        !> \brief Exponential PPFD regression coefficient C for sunlit leaves at level 5
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 5 (bottom)
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_5_SUN     =  1.090_rk

        !> \brief Exponential PPFD regression coefficient D for sunlit leaves at level 1
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 1
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_1_SUN     =  0.002_rk

        !> \brief Exponential PPFD regression coefficient D for sunlit leaves at level 2
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 2
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_2_SUN     =  -0.128_rk

        !> \brief Exponential PPFD regression coefficient D for sunlit leaves at level 3
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 3
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_3_SUN     =  -0.298_rk

        !> \brief Exponential PPFD regression coefficient D for sunlit leaves at level 4
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 4
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_4_SUN     =  -0.445_rk

        !> \brief Exponential PPFD regression coefficient D for sunlit leaves at level 5
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 5
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_5_SUN     =  -0.535_rk

        !> \brief Exponential PPFD regression coefficient C for shaded leaves at level 1
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 1
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_1_SHADE   =  0.871_rk

        !> \brief Exponential PPFD regression coefficient C for shaded leaves at level 2
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 2
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_2_SHADE   =  0.890_rk

        !> \brief Exponential PPFD regression coefficient C for shaded leaves at level 3
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 3
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_3_SHADE   =  0.916_rk

        !> \brief Exponential PPFD regression coefficient C for shaded leaves at level 4
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 4
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_4_SHADE   =  0.941_rk

        !> \brief Exponential PPFD regression coefficient C for shaded leaves at level 5
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 5
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: CTEMP_5_SHADE   =  0.956_rk

        !> \brief Exponential PPFD regression coefficient D for shaded leaves at level 1
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 1
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_1_SHADE   =  0.015_rk

        !> \brief Exponential PPFD regression coefficient D for shaded leaves at level 2
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 2
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_2_SHADE   =  -0.141_rk

        !> \brief Exponential PPFD regression coefficient D for shaded leaves at level 3
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 3
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_3_SHADE   =  -0.368_rk
        !> \brief Exponential PPFD regression coefficient D for shaded leaves at level 4
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 4
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_4_SHADE   =  -0.592_rk

        !> \brief Exponential PPFD regression coefficient D for shaded leaves at level 5
        !> \details Regression coefficient from Silva et al. (2020) Table 1, level 5
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: DTEMP_5_SHADE   =  -0.743_rk

        !> \brief Fraction of PAR in solar irradiance
        !> \details Fraction of incoming solar irradiance that is photosynthetically active radiation
        !! \param units dimensionless
        REAL(RK),          PARAMETER     :: FRAC_PAR        =  0.5_rk

        !> \brief Regression coefficient C for sun leaves
        !> \details Height-interpolated regression coefficient C for sun leaves
        REAL(RK) :: CTEMP_SUN(SIZE(ZK))

        !> \brief Regression coefficient D for sun leaves
        !> \details Height-interpolated regression coefficient D for sun leaves
        REAL(RK) :: DTEMP_SUN(SIZE(ZK))

        !> \brief Regression coefficient C for shade leaves
        !> \details Height-interpolated regression coefficient C for shade leaves
        REAL(RK) :: CTEMP_SHADE(SIZE(ZK))

        !> \brief Regression coefficient D for shade leaves
        !> \details Height-interpolated regression coefficient D for shade leaves
        REAL(RK) :: DTEMP_SHADE(SIZE(ZK))

        !> \brief Loop index
        !> \details Loop index for height levels
        integer i

! Use exponential PPFD model based on Silva et al. (2020) to get approx. sun/shade PPFD
! through canopy
!Citation:
!Silva, S. J., Heald, C. L., and Guenther, A. B.: Development of a reduced-complexity plant canopy
!physics surrogate model for use in chemical transport models: a case study with GEOS-Chem v12.3.0,
!Geosci. Model Dev., 13, 2569–2585, https://doi.org/10.5194/gmd-13-2569-2020, 2020.
        do i=1, SIZE(ZK)  !calculate linear change in parameters interpolated to Silva et al. 5 layer canopy regions
            if (ZK(i) .gt. FCH) then ! above canopy, PPFD_leaf = PPFD_toc (toc=top of canopy)
                CTEMP_SUN(i)   = 0.0
                DTEMP_SUN(i)   = 0.0
                CTEMP_SHADE(i) = 0.0
                DTEMP_SHADE(i) = 0.0
            else if (ZK(i) .le. FCH .and. ZK(i) .gt. FCH*(4.0_rk/5.0_rk)) then  !Level 1 - 2
                CTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ CTEMP_2_SUN,CTEMP_1_SUN /),ZK(i))
                DTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ DTEMP_2_SUN,DTEMP_1_SUN /),ZK(i))
                CTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ CTEMP_2_SHADE,CTEMP_1_SHADE /),ZK(i))
                DTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ DTEMP_2_SHADE,DTEMP_1_SHADE /),ZK(i))
            else if (ZK(i) .le. FCH*(4.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(3.0_rk/5.0_rk)) then  !Level 2 - 3
                CTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ CTEMP_3_SUN,CTEMP_2_SUN /),ZK(i))
                DTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ DTEMP_3_SUN,DTEMP_2_SUN /),ZK(i))
                CTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ CTEMP_3_SHADE,CTEMP_2_SHADE /),ZK(i))
                DTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ DTEMP_3_SHADE,DTEMP_2_SHADE /),ZK(i))
            else if (ZK(i) .le. FCH*(3.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(2.0_rk/5.0_rk)) then  !Level 3 - 4
                CTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ CTEMP_4_SUN,CTEMP_3_SUN /),ZK(i))
                DTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ DTEMP_4_SUN,DTEMP_3_SUN /),ZK(i))
                CTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ CTEMP_4_SHADE,CTEMP_3_SHADE /),ZK(i))
                DTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ DTEMP_4_SHADE,DTEMP_3_SHADE /),ZK(i))
            else if (ZK(i) .le. FCH*(2.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                CTEMP_SUN(i)   = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ CTEMP_5_SUN,CTEMP_4_SUN /),ZK(i))
                DTEMP_SUN(i)   = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ DTEMP_5_SUN,DTEMP_4_SUN /),ZK(i))
                CTEMP_SHADE(i) = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ CTEMP_5_SHADE,CTEMP_4_SHADE /),ZK(i))
                DTEMP_SHADE(i) = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ DTEMP_5_SHADE,DTEMP_4_SHADE /),ZK(i))
            end if
        end do

        PPFD_SUN     = FRAC_PAR * SFCRAD * EXP(CTEMP_SUN + DTEMP_SUN * LAI)  !Silva et al. W/m2 --> umol m-2 s-1
        PPFD_SHADE   = FRAC_PAR * SFCRAD * EXP(CTEMP_SHADE + DTEMP_SHADE * LAI)
        PPFD_AVE = (PPFD_SUN*FSUN) + (PPFD_SHADE*(1.0-FSUN)) ! average = sum sun and shade weighted by sunlit fraction

    END SUBROUTINE CANOPY_PPFD_EXP

end module canopy_rad_mod
