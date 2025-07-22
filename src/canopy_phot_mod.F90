!> \file canopy_phot_mod.F90
!! \brief Canopy Photolysis Module
!! \details This module contains subroutines for calculating photolysis attenuation
!! within forest canopies based on the algorithms described in Makar et al. (2017).
!! The module computes how photolysis rates are reduced within the canopy due to
!! shading by leaves and branches.
!!
!! \author Patrick C. Campbell
!! \date June 2022
!!
!! \references
!! Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy
!! shading and turbulence on boundary layer ozone. Nat Commun 8, 15243 (2017).
!! https://doi.org/10.1038/ncomms15243

!> \defgroup phot_mod Canopy Photolysis Module
!! \brief Module for photolysis attenuation calculations in canopies
!! \{

module canopy_phot_mod

    implicit none

contains

!> \brief Calculate photolysis attenuation in canopy
!! \details This subroutine computes the photolysis correction factors within a forest canopy
!! using the exponential attenuation model from Makar et al. (2017). The calculation
!! accounts for:
!! - Fractional cumulative leaf area index (FCLAI) profiles
!! - Total leaf area index (LAI)
!! - Clumping index to account for non-random leaf distribution
!! - Solar zenith angle effects on light penetration
!!
!! The photolysis correction factor represents the fraction of photolysis rates
!! relative to above-canopy conditions at each height within the canopy.
!!
!! \param[in] FCLAI Fractional cumulative LAI shapes of plant surface distribution (nondimensional)
!! \param[in] LAI Model input total Leaf Area Index (m²/m²)
!! \param[in] CLU Model input Clumping Index (nondimensional)
!! \param[in] COSZEN Model input Cosine Solar Zenith Angle (nondimensional)
!! \param[out] RJCF Photolysis correction factor (nondimensional)
    SUBROUTINE CANOPY_PHOT( FCLAI, LAI, CLU, COSZEN, RJCF )
        use canopy_const_mod, ONLY: rk     !> constants for canopy models

!> \defgroup phot_inputs Input Variables
!! \brief Input parameters for photolysis calculations
!! \{
        REAL(RK),    INTENT( IN )  :: FCLAI(:)           !> Model input Fractional (z) shapes of the
        !! plant surface distribution (nondimensional), i.e., a Fractional Cumulative LAI
        REAL(RK),    INTENT( IN )  :: LAI             !> Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )  :: CLU             !> Model input Clumping Index
        REAL(RK),    INTENT( IN )  :: COSZEN          !> Model input Cosine Solar Zenith Angle
!> \}

!> \defgroup phot_outputs Output Variables
!! \brief Output parameters for photolysis calculations
!! \{
        REAL(RK),    INTENT( OUT ) :: RJCF(:)          !> Photolysis correction factor
!> \}

!> \brief Calculate photolysis correction factor using exponential attenuation
!! \details Uses Eq. 1 from Makar et al. (2017) to compute the exponential attenuation
!! of photolysis rates through the canopy based on cumulative leaf area and solar angle.
        RJCF = MAX(1.0E-10_rk, EXP(-1.0_rk*(0.5_rk*(LAI*(1.0_rk-FCLAI))*CLU)/MAX(0.05_rk, COSZEN)))

    END SUBROUTINE CANOPY_PHOT

!> \}

end module canopy_phot_mod
