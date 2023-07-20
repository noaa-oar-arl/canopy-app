module canopy_rad_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

end module canopy_rad_mod
