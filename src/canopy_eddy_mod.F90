!> \file canopy_eddy_mod.F90
!> \brief Canopy eddy diffusivity calculations
!> \details This module contains routines for calculating eddy diffusivity (Kz)
!>          within and just above forest canopies using parameterizations
!>          based on Makar et al. (2017) and Raupach (1989).
!> \author P. C. Campbell
!> \date Jun 2022

!> \defgroup canopy_eddy Canopy Eddy Diffusivity
!> \brief Eddy diffusivity calculations within and above canopy

module canopy_eddy_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !> \brief Calculate canopy eddy diffusivity
    !> \details Computes Eddy Diffusivity, Kz, within and just above canopy using
    !!          stability-dependent parameterizations for different atmospheric conditions
    !> \ingroup canopy_eddy
    !> \param[in] HCM Height of canopy top (m)
    !> \param[in] ZK Above/Below canopy height, z (m)
    !> \param[in] USTAR Model input friction velocity (m/s)
    !> \param[in] MOL Model input Monin-Obukhov Length (m)
    !> \param[out] KZ Estimated Eddy Diffusivity with canopy turbulence (m²/s)
    !> \author P. C. Campbell
    !> \date Jun 2022
    !> \note Based on Makar et al. (2017) algorithms for canopy turbulence
    !> \cite Makar, P., Staebler, R., Akingunola, A. et al. (2017).
    !!       The effects of forest canopy shading and turbulence on boundary layer ozone.
    !!       Nature Communications, 8, 15243. https://doi.org/10.1038/ncomms15243
    !> \cite Raupach, M. R. (1989). A Practical Lagrangian method for relating scalar
    !!       concentrations to source distributions in vegetation canopies.
    !!       Quarterly Journal of the Royal Meteorological Society, 115, 609-632.
    SUBROUTINE CANOPY_EDDYX( HCM, ZK, USTAR, MOL, KZ )

!-----------------------------------------------------------------------

! Description:
!     computes Eddy Diffusivity, Kz, within and just above canopy.

! Preconditions:
!     in-canopy height, and model friction velocity and MOL

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Makar et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone vertical diffusion calculation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk, pi !constants for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: ZK              ! Above/Below canopy height, z (m)
        REAL(RK),    INTENT( IN )  :: USTAR           ! Model input friction velocity (m/s)
        REAL(RK),    INTENT( IN )  :: MOL             ! Model input Monin-Obukhov Length
        REAL(RK),    INTENT( OUT ) :: KZ              ! Estimated Eddy Diffusivity with canopy turbulence

!     Local variables
        !> \brief Local canopy stability parameter
        !> \details Canopy stability parameter (hc/MOL)
        !! \param units dimensionless
        real(rk)                   :: hol

        !> \brief Turbulence length scale in canopy
        !> \details Turbulence length scale in canopy from Raupach (1989)
        !! \param units meters (m)
        real(rk)                   :: tlc

        !> \brief Sigma parameter R
        !> \details Sigma parameter R from Makar et al. (2017)
        !! \param units dimensionless
        real(rk)                   :: rr

        !> \brief Sigma parameter A
        !> \details Sigma parameter A from Makar et al. (2017)
        !! \param units dimensionless
        real(rk)                   :: aa

        !> \brief Sigma parameter B
        !> \details Sigma parameter B from Makar et al. (2017)
        !! \param units dimensionless
        real(rk)                   :: bb

        !> \brief Eulerian vertical velocity
        !> \details Eulerian vertical velocity variance
        !! \param units m/s
        real(rk)                   :: sigma

! Citation:
!Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone.
!Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243
!Raupauch M. R. A Practical Lagrangian method for relating scalar concentrations to
! source distributions in vegetation canopies. Q. J. R. Meteor. Soc. (1989), 115, pp 609-632
! Eqs. 2-9 on pgs 10-11.

        hol = HCM/MOL
        tlc = (HCM/USTAR) * (                        &
            (0.256_rk * (ZK-(0.75_rk*HCM))/HCM ) +      &
            (0.492_rk*EXP((-0.256_rk*ZK/HCM)/0.492_rk)) )
        sigma = 0.0_rk

        IF ( hol .LT. -0.1 )  THEN  !UNSTABLE
            IF ( ZK/HCM .GT. 1.25_rk ) THEN !SIGMACAN = Eulerian vertical velocity variance
                sigma = 1.25_rk*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                sigma = USTAR * ( 0.75_rk + (0.5_rk * COS((pi/1.06818_rk) *     &
                    (1.25_rk - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25_rk*USTAR
            END IF
        END IF

        IF ( hol .GE. -0.1  .AND. hol .LT. 0.1 )   THEN  !NEUTRAL
            IF ( ZK/HCM .GT. 1.25 ) THEN
                sigma = 1.0_rk*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                sigma = USTAR * ( 0.625_rk + (0.375_rk * COS((pi/1.06818_rk) *  &
                    (1.25_rk - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25_rk*USTAR
            END IF
        END IF

        IF ( hol .GE.  0.1  .AND.  hol .LT. 0.9 )   THEN  !STABLE
            IF ( ZK/HCM .GT. 1.25 ) THEN
                sigma = 0.25_rk*(4.375_rk - (3.75_rk*hol))*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                rr=4.375_rk-(3.75_rk*hol)
                aa=(0.125_rk*rr) + 0.125_rk
                bb=(0.125_rk*rr) - 0.125_rk
                sigma = USTAR * ( aa + (bb * COS((pi/1.06818_rk) *   &
                    (1.25_rk - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25_rk*USTAR
            END IF
        END IF

        IF ( hol .GE.  0.9 ) THEN  !VERY STABLE
            sigma = 0.25_rk*USTAR
        END IF

        KZ = (sigma*sigma)*tlc

    END SUBROUTINE CANOPY_EDDYX

end module canopy_eddy_mod
