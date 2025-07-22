!> \file canopy_tleaf_mod.F90
!> \brief Leaf temperature calculation module for canopy model
!> \details This module implements linear interpolation methods for calculating
!!          leaf temperatures for sunlit and shaded leaves throughout the canopy
!!          based on Silva et al. (2020) algorithms.
!> \author P.C. Campbell
!> \date June 2023
!> \version 1.0

!> \defgroup tleaf_group Leaf Temperature Calculations
!> \brief Routines for calculating leaf temperatures in forest canopies
!> \details This group contains subroutines for calculating leaf temperatures
!!          for sunlit and shaded leaves using linear interpolation methods
!!          based on Silva et al. (2020) parameterizations.
!> \{

module canopy_tleaf_mod

    implicit none

contains

!> \brief Linear interpolation method for leaf temperature calculation in canopy
!> \details Computes linear interpolation method for leaf temperature (sun/shade)
!!          in canopy based on Silva et al. (2020) algorithms. Uses linearized
!!          2-m temperature to leaf temperature parameters with regression
!!          coefficients for 5 canopy layers.
!> \param ZK Input model heights (m)
!> \param FCH Model input canopy height (m)
!> \param TEMP2 Model input 2-m temperature (K)
!> \param FSUN Sunlit/shaded fraction from photolysis correction factor
!> \param TLEAF_SUN Leaf temperature for sunlit leaves (K) [output]
!> \param TLEAF_SHADE Leaf temperature for shaded leaves (K) [output]
!> \param TLEAF_AVE Average leaf temperature for sun/shaded leaves (K) [output]
!> \author P.C. Campbell
!> \date June 2023
!> \note Based on Silva et al. (2020) - Development of a reduced-complexity
!!       plant canopy physics surrogate model for use in chemical transport models
    SUBROUTINE CANOPY_TLEAF_LIN( ZK, FCH, TEMP2, FSUN, &
        TLEAF_SUN, TLEAF_SHADE, TLEAF_AVE)

        use canopy_const_mod, ONLY: RK     !< Constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal !< Linear interpolation utility

        !> \name Input Parameters
        !> \{
        REAL(RK),    INTENT( IN )       :: ZK(:)                          !< Input model heights (m)
        REAL(RK),    INTENT( IN )       :: FCH                            !< Model input canopy height (m)
        REAL(RK),    INTENT( IN )       :: TEMP2                          !< Model input 2-m Temperature (K)
        REAL(RK),    INTENT( IN )       :: FSUN(:)                        !< Sunlit/Shaded fraction from photolysis correction factor
        !> \}

        !> \name Output Parameters
        !> \{
        REAL(RK),    INTENT( OUT )      :: TLEAF_SUN(SIZE(ZK))            !< Leaf temp for sunlit leaves (K)
        REAL(RK),    INTENT( OUT )      :: TLEAF_SHADE(SIZE(ZK))          !< Leaf temp for shaded leaves (K)
        REAL(RK),    INTENT( OUT )      :: TLEAF_AVE(SIZE(ZK))            !< Ave Leaf temp for sun/shaded leaves (K)
        !> \}

        !> \name Local Variables - Silva et al. (2020) Regression Coefficients
        !> \details Linearized 2-m temp --> leaf temp parameters based on Table 1 in Silva et al. (2020)
        !> \{
        !> Regression coefficient A for sunlit leaves at canopy level 1 (top)
        REAL(RK),          PARAMETER     :: ATEMP_1_SUN     =  -13.891_rk
        !> Regression coefficient A for sunlit leaves at canopy level 2
        REAL(RK),          PARAMETER     :: ATEMP_2_SUN     =  -12.322_rk
        !> Regression coefficient A for sunlit leaves at canopy level 3
        REAL(RK),          PARAMETER     :: ATEMP_3_SUN     =  -1.032_rk
        !> Regression coefficient A for sunlit leaves at canopy level 4
        REAL(RK),          PARAMETER     :: ATEMP_4_SUN     =  -5.172_rk
        !> Regression coefficient A for sunlit leaves at canopy level 5 (bottom)
        REAL(RK),          PARAMETER     :: ATEMP_5_SUN     =  -5.589_rk
        !> Regression coefficient B for sunlit leaves at canopy level 1 (top)
        REAL(RK),          PARAMETER     :: BTEMP_1_SUN     =  1.064_rk
        !> Regression coefficient B for sunlit leaves at canopy level 2
        REAL(RK),          PARAMETER     :: BTEMP_2_SUN     =  1.057_rk
        !> Regression coefficient B for sunlit leaves at canopy level 3
        REAL(RK),          PARAMETER     :: BTEMP_3_SUN     =  1.031_rk
        !> Regression coefficient B for sunlit leaves at canopy level 4
        REAL(RK),          PARAMETER     :: BTEMP_4_SUN     =  1.050_rk
        !> Regression coefficient B for sunlit leaves at canopy level 5 (bottom)
        REAL(RK),          PARAMETER     :: BTEMP_5_SUN     =  1.051_rk
        !> Regression coefficient A for shaded leaves at canopy level 1 (top)
        REAL(RK),          PARAMETER     :: ATEMP_1_SHADE   =  -12.846_rk
        !> Regression coefficient A for shaded leaves at canopy level 2
        REAL(RK),          PARAMETER     :: ATEMP_2_SHADE   =  -11.343_rk
        !> Regression coefficient A for shaded leaves at canopy level 3
        REAL(RK),          PARAMETER     :: ATEMP_3_SHADE   =  -1.068_rk
        !> Regression coefficient A for shaded leaves at canopy level 4
        REAL(RK),          PARAMETER     :: ATEMP_4_SHADE   =  -5.551_rk
        !> Regression coefficient A for shaded leaves at canopy level 5 (bottom)
        REAL(RK),          PARAMETER     :: ATEMP_5_SHADE   =  -5.955_rk
        !> Regression coefficient B for shaded leaves at canopy level 1 (top)
        REAL(RK),          PARAMETER     :: BTEMP_1_SHADE   =  1.060_rk
        !> Regression coefficient B for shaded leaves at canopy level 2
        REAL(RK),          PARAMETER     :: BTEMP_2_SHADE   =  1.053_rk
        !> Regression coefficient B for shaded leaves at canopy level 3
        REAL(RK),          PARAMETER     :: BTEMP_3_SHADE   =  1.031_rk
        !> Regression coefficient B for shaded leaves at canopy level 4
        REAL(RK),          PARAMETER     :: BTEMP_4_SHADE   =  1.051_rk
        !> Regression coefficient B for shaded leaves at canopy level 5 (bottom)
        REAL(RK),          PARAMETER     :: BTEMP_5_SHADE   =  1.053_rk
        !> \}

        !> \name Working Arrays
        !> \{
        REAL(RK) :: ATEMP_SUN(SIZE(ZK))            !< Regression coefficient A for sun leaves (Silva et al., 2020)
        REAL(RK) :: BTEMP_SUN(SIZE(ZK))            !< Regression coefficient B for sun leaves
        REAL(RK) :: ATEMP_SHADE(SIZE(ZK))          !< Regression coefficient A for shade leaves
        REAL(RK) :: BTEMP_SHADE(SIZE(ZK))          !< Regression coefficient B for shade leaves
        !> \}

        integer i !< Loop index

        !> Use linear canopy temperature model based on Silva et al. (2020) to get approx. sun/shade leaf temperatures
        !! through canopy (ignores effect of wind speed on leaf boundary layer ~ 1 % error/bias)
        !! Citation:
        !! Silva, S. J., Heald, C. L., and Guenther, A. B.: Development of a reduced-complexity plant canopy
        !! physics surrogate model for use in chemical transport models: a case study with GEOS-Chem v12.3.0,
        !! Geosci. Model Dev., 13, 2569–2585, https://doi.org/10.5194/gmd-13-2569-2020, 2020.
        do i=1, SIZE(ZK)  !< Calculate linear change in parameters interpolated to Silva et al. 5 layer canopy regions
            if (ZK(i) .gt. FCH) then !< Above canopy, Tleaf = Tair
                ATEMP_SUN(i)   = 0.0
                BTEMP_SUN(i)   = 1.0
                ATEMP_SHADE(i) = 0.0
                BTEMP_SHADE(i) = 1.0
            else if (ZK(i) .le. FCH .and. ZK(i) .gt. FCH*(4.0_rk/5.0_rk)) then  !< Level 1 - 2
                ATEMP_SUN(i)   = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ ATEMP_2_SUN,ATEMP_1_SUN /),ZK(i))
                BTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ BTEMP_2_SUN,BTEMP_1_SUN /),ZK(i))
                ATEMP_SHADE(i) = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ ATEMP_2_SHADE,ATEMP_1_SHADE /),ZK(i))
                BTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(4.0_rk/5.0_rk),FCH /), &
                    (/ BTEMP_2_SHADE,BTEMP_1_SHADE /),ZK(i))
            else if (ZK(i) .le. FCH*(4.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(3.0_rk/5.0_rk)) then  !< Level 2 - 3
                ATEMP_SUN(i)   = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ ATEMP_3_SUN,ATEMP_2_SUN /),ZK(i))
                BTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ BTEMP_3_SUN,BTEMP_2_SUN /),ZK(i))
                ATEMP_SHADE(i) = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ ATEMP_3_SHADE,ATEMP_2_SHADE /),ZK(i))
                BTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(3.0_rk/5.0_rk),FCH*(4.0_rk/5.0_rk) /), &
                    (/ BTEMP_3_SHADE,BTEMP_2_SHADE /),ZK(i))
            else if (ZK(i) .le. FCH*(3.0_rk/5.0_rk) .and. ZK(i) .gt. FCH*(2.0_rk/5.0_rk)) then  !< Level 3 - 4
                ATEMP_SUN(i)   = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ ATEMP_4_SUN,ATEMP_3_SUN /),ZK(i))
                BTEMP_SUN(i)   = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ BTEMP_4_SUN,BTEMP_3_SUN /),ZK(i))
                ATEMP_SHADE(i) = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ ATEMP_4_SHADE,ATEMP_3_SHADE /),ZK(i))
                BTEMP_SHADE(i) = interp_linear1_internal((/ FCH*(2.0_rk/5.0_rk),FCH*(3.0_rk/5.0_rk) /), &
                    (/ BTEMP_4_SHADE,BTEMP_3_SHADE /),ZK(i))
            else if (ZK(i) .le. FCH*(2.0_rk/5.0_rk) ) then  !< Level 4 - Bottom
                ATEMP_SUN(i)   = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ ATEMP_5_SUN,ATEMP_4_SUN /),ZK(i))
                BTEMP_SUN(i)   = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ BTEMP_5_SUN,BTEMP_4_SUN /),ZK(i))
                ATEMP_SHADE(i) = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ ATEMP_5_SHADE,ATEMP_4_SHADE /),ZK(i))
                BTEMP_SHADE(i) = interp_linear1_internal((/ ZK(1),FCH*(2.0_rk/5.0_rk) /), &
                    (/ BTEMP_5_SHADE,BTEMP_4_SHADE /),ZK(i))
            end if
        end do

        !> Calculate final leaf temperatures using regression coefficients
        TLEAF_SUN   = ATEMP_SUN + (BTEMP_SUN*TEMP2)
        TLEAF_SHADE = ATEMP_SHADE + (BTEMP_SHADE*TEMP2)
        !> Average = sum sun and shade weighted by sunlit fraction
        TLEAF_AVE = (TLEAF_SUN*FSUN) + (TLEAF_SHADE*(1.0-FSUN))

    END SUBROUTINE CANOPY_TLEAF_LIN

!> \}

end module canopy_tleaf_mod
