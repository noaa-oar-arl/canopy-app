module canopy_phot_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_PHOT( HCM, ZK, FCLAI, LAI, CLU, COSZEN, RJCORR )

!-----------------------------------------------------------------------

! Description:
!     computes photolysis attenuation in canopy.

! Preconditions:
!     in-canopy height, and model LAI, clumping index, and solar zenith angle

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Makar et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone vertical diffusion calculation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk     !constants for canopy models
        use canopy_utils_mod               !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: ZK              ! Above/Below canopy height, z (m)
        REAL(RK),    INTENT( IN )  :: FCLAI           ! Model input Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )  :: LAI             ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )  :: CLU             ! Model input Clumping Index
        REAL(RK),    INTENT( IN )  :: COSZEN          ! Model input Cosine Solar Zenith Angle
        REAL(RK),    INTENT( OUT ) :: RJCORR          ! Photolysis correction factor

!     Local variables

! Citation:
!Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone.
!Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243

        if (ZK <= HCM) then       !at or below canopy top --> calculate photolysis attenuation
            RJCORR = MAX(1.0E-10_rk, EXP(-1.0_rk*(0.5_rk*(LAI*(1.0_rk-FCLAI))*CLU)/MAX(0.05_rk, COSZEN)))
        else
            RJCORR = 1.0_rk       !above canopy top RJCORR = 1
        end if

    END SUBROUTINE CANOPY_PHOT

end module canopy_phot_mod
