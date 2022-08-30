module canopy_eddyx_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
        use canopy_utils_mod               !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: ZK              ! Above/Below canopy height, z (m)
        REAL(RK),    INTENT( IN )  :: USTAR           ! Model input friction velocity (m/s)
        REAL(RK),    INTENT( IN )  :: MOL             ! Model input Monin-Obukhov Length
        REAL(RK),    INTENT( OUT ) :: KZ              ! Estimated Eddy Diffusivity with canopy turbulence

!     Local variables
        real(rk)                   :: hol             ! local canopy stability parameter (hc/MOL)
        real(rk)                   :: tlc             ! turbulence length scale in canopy (m)
        real(rk)                   :: rr              ! sigma parameter R (Makar et al., 2017)
        real(rk)                   :: aa              ! sigma parameter A (Makar et al., 2017)
        real(rk)                   :: bb              ! sigma parameter B (Makar et al., 2017)
        real(rk)                   :: sigma           ! eulerian vertial velocity (m/s)

! Citation:
!Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone.
!Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243
!Raupauch M. R. A Practical Lagrangian method for relating scalar concentrations to
! source distributions in vegetation canopies. Q. J. R. Meteor. Soc. (1989), 115, pp 609-632

        hol = HCM/MOL
        tlc = (HCM/USTAR) * (                        &
            (0.256 * (ZK-(0.75*HCM))/HCM ) +      &
            (0.492*EXP((-0.256*ZK/HCM)/0.492)) )
        sigma = 0.0_rk

        IF ( hol .LT. -0.1 )  THEN  !UNSTABLE
            IF ( ZK/HCM .GT. 1.25 ) THEN !SIGMACAN = Eulerian vertical velocity variance
                sigma = 1.25*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                sigma = USTAR * ( 0.75 + (0.5 * COS((pi/1.06818) *     &
                    (1.25 - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25*USTAR
            END IF
        END IF

        IF ( hol .GE. -0.1  .AND. hol .LT. 0.1 )   THEN  !NEUTRAL
            IF ( ZK/HCM .GT. 1.25 ) THEN
                sigma = 1.0*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                sigma = USTAR * ( 0.625 + (0.375 * COS((pi/1.06818) *  &
                    (1.25 - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25*USTAR
            END IF
        END IF

        IF ( hol .GE.  0.1  .AND.  hol .LT. 0.9 )   THEN  !STABLE
            IF ( ZK/HCM .GT. 1.25 ) THEN
                sigma = 0.25*(4.375 - (3.75*hol))*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                rr=4.375-(3.75*hol)
                aa=(0.125*rr) + 0.125
                bb=(0.125*rr) - 0.125
                sigma = USTAR * ( aa + (bb * COS((pi/1.06818) *   &
                    (1.25 - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25*USTAR
            END IF
        END IF

        IF ( hol .GE.  0.9 ) THEN  !VERY STABLE
            sigma = 0.25*USTAR
        END IF

        KZ = (sigma*sigma)*tlc

!       print*, 'ZK/HCM = ', ZK/HCM, 'KZ = ', KZ




    END SUBROUTINE CANOPY_EDDYX

end module canopy_eddyx_mod
