

# File canopy\_eddy\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_eddy\_mod.F90**](canopy__eddy__mod_8F90.md)

[Go to the documentation of this file](canopy__eddy__mod_8F90.md)


```Fortran



module canopy_eddy_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_eddyx( HCM, ZK, USTAR, MOL, KZ )

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
        real(rk)                   :: hol

        real(rk)                   :: tlc

        real(rk)                   :: rr

        real(rk)                   :: aa

        real(rk)                   :: bb

        real(rk)                   :: sigma

! Citation:
!Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone.
!Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243
!Raupauch M. R. A Practical Lagrangian method for relating scalar concentrations to
! source distributions in vegetation canopies. Q. J. R. Meteor. Soc. (1989), 115, pp 609-632
! Eqs. 2-9 on pgs 10-11.

        hol = hcm/mol
        tlc = (hcm/ustar) * (                        &
            (0.256_rk * (zk-(0.75_rk*hcm))/hcm ) +      &
            (0.492_rk*exp((-0.256_rk*zk/hcm)/0.492_rk)) )
        sigma = 0.0_rk

        IF ( hol .LT. -0.1 )  THEN  !UNSTABLE
            IF ( zk/hcm .GT. 1.25_rk ) THEN !SIGMACAN = Eulerian vertical velocity variance
                sigma = 1.25_rk*ustar
            END IF
            IF ( zk/hcm .GE. 0.175  .AND.  zk/hcm .LE. 1.25 ) THEN
                sigma = ustar * ( 0.75_rk + (0.5_rk * cos((pi/1.06818_rk) *     &
                    (1.25_rk - (zk/hcm)))) )
            END IF
            IF ( zk/hcm .LT. 0.175 )  THEN
                sigma = 0.25_rk*ustar
            END IF
        END IF

        IF ( hol .GE. -0.1  .AND. hol .LT. 0.1 )   THEN  !NEUTRAL
            IF ( zk/hcm .GT. 1.25 ) THEN
                sigma = 1.0_rk*ustar
            END IF
            IF ( zk/hcm .GE. 0.175  .AND.  zk/hcm .LE. 1.25 ) THEN
                sigma = ustar * ( 0.625_rk + (0.375_rk * cos((pi/1.06818_rk) *  &
                    (1.25_rk - (zk/hcm)))) )
            END IF
            IF ( zk/hcm .LT. 0.175 )  THEN
                sigma = 0.25_rk*ustar
            END IF
        END IF

        IF ( hol .GE.  0.1  .AND.  hol .LT. 0.9 )   THEN  !STABLE
            IF ( zk/hcm .GT. 1.25 ) THEN
                sigma = 0.25_rk*(4.375_rk - (3.75_rk*hol))*ustar
            END IF
            IF ( zk/hcm .GE. 0.175  .AND.  zk/hcm .LE. 1.25 ) THEN
                rr=4.375_rk-(3.75_rk*hol)
                aa=(0.125_rk*rr) + 0.125_rk
                bb=(0.125_rk*rr) - 0.125_rk
                sigma = ustar * ( aa + (bb * cos((pi/1.06818_rk) *   &
                    (1.25_rk - (zk/hcm)))) )
            END IF
            IF ( zk/hcm .LT. 0.175 )  THEN
                sigma = 0.25_rk*ustar
            END IF
        END IF

        IF ( hol .GE.  0.9 ) THEN  !VERY STABLE
            sigma = 0.25_rk*ustar
        END IF

        kz = (sigma*sigma)*tlc

    END SUBROUTINE canopy_eddyx

end module canopy_eddy_mod
```


