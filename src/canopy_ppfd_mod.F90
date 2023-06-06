module canopy_ppfd_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_PPFD_EXP( ZK, FCH, SFCRAD, LAI, PPFD_SUN, PPFD_SHADE)

!-----------------------------------------------------------------------

! Description:
!     computes linear interpolation method for tleaf sun/shade in canopy.

! Preconditions:
!     in-canopy height, and model LAI, clumping index, and solar zenith angle

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/23 by PCC
!     Jun 2023 P.C. Campbell: Initial standalone tleaf linear subroutine based on
!                             Silva et al. (2020) algorithms
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
        REAL(RK),    INTENT( OUT )      :: PPFD_SUN(SIZE(ZK))             ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( OUT )      :: PPFD_SHADE(SIZE(ZK))           ! PPFD for shaded leaves (umol phot/m2 s)

!      LOCAL
        REAL(RK),          PARAMETER     :: CTEMP_1_SUN     =  1.083_rk   !Exponential 2-m PPFD --> PPFD parameters (Level 1 =
        !top of canopy
        REAL(RK),          PARAMETER     :: CTEMP_2_SUN     =  1.096_rk   !Based on Table 1 in Silva et al. (2022)
        REAL(RK),          PARAMETER     :: CTEMP_3_SUN     =  1.104_rk   !
        REAL(RK),          PARAMETER     :: CTEMP_4_SUN     =  1.098_rk   !
        REAL(RK),          PARAMETER     :: CTEMP_5_SUN     =  1.090_rk   !...
        REAL(RK),          PARAMETER     :: DTEMP_1_SUN     =  0.002_rk   !...
        REAL(RK),          PARAMETER     :: DTEMP_2_SUN     =  -0.128_rk  !...
        REAL(RK),          PARAMETER     :: DTEMP_3_SUN     =  -0.298_rk  !...
        REAL(RK),          PARAMETER     :: DTEMP_4_SUN     =  -0.445_rk  !...
        REAL(RK),          PARAMETER     :: DTEMP_5_SUN     =  -0.535_rk  !...
        REAL(RK),          PARAMETER     :: CTEMP_1_SHADE   =  0.871_rk   !...
        REAL(RK),          PARAMETER     :: CTEMP_2_SHADE   =  0.890_rk   !...
        REAL(RK),          PARAMETER     :: CTEMP_3_SHADE   =  0.916_rk   !...
        REAL(RK),          PARAMETER     :: CTEMP_4_SHADE   =  0.941_rk   !...
        REAL(RK),          PARAMETER     :: CTEMP_5_SHADE   =  0.956_rk   !...
        REAL(RK),          PARAMETER     :: DTEMP_1_SHADE   =  0.015_rk   !...
        REAL(RK),          PARAMETER     :: DTEMP_2_SHADE   =  -0.141_rk  !...
        REAL(RK),          PARAMETER     :: DTEMP_3_SHADE   =  -0.368_rk  !...
        REAL(RK),          PARAMETER     :: DTEMP_4_SHADE   =  -0.592_rk  !...
        REAL(RK),          PARAMETER     :: DTEMP_5_SHADE   =  -0.743_rk  !...

        REAL(RK),          PARAMETER     :: FRAC_PAR        =  0.5_rk     !Fraction of incoming solar irradiance that is PAR

        REAL(RK) :: CTEMP_SUN(SIZE(ZK))            ! Regression coefficient C for sun leaves
        REAL(RK) :: DTEMP_SUN(SIZE(ZK))            ! Regression coefficient D for sun leaves
        REAL(RK) :: CTEMP_SHADE(SIZE(ZK))          ! Regression coefficient C for shade leaves
        REAL(RK) :: DTEMP_SHADE(SIZE(ZK))          ! Regression coefficient D for shade leaves

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

    END SUBROUTINE CANOPY_PPFD_EXP

end module canopy_ppfd_mod
