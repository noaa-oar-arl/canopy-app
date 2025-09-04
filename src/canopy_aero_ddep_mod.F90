!> \file canopy_aero_ddep_mod.F90
!> \brief Aerosol dry deposition calculations for canopy model (sub-canopy)
!> \details Implements sub-canopy aerosol dry deposition following Katul et al. (2010)
!!          https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JD012853
!!          and issue #148. Calculates deposition velocities for aerosols in multilayer canopy.
!> \author NOAA ARL
!> \date September 2025
!> \version 1.0

module canopy_aero_ddep_mod

    implicit none
    private
    public :: canopy_aero_ddep_katul2010

contains

!> \brief Sub-canopy aerosol dry deposition velocity following Katul et al. (2010)
!> \param nlev Number of canopy layers
!> \param z Array of heights above ground (m)
!> \param h_can Canopy height (m)
!> \param lad Leaf area density profile (m^2/m^3)
!> \param u Array of wind speed profile (m/s)
!> \param d_p Aerosol particle diameter (m)
!> \param rho_p Particle density (kg/m^3)
!> \param T Air temperature profile (K)
!> \param P Air pressure profile (Pa)
!> \param vdep_aero Output: aerosol deposition velocity profile (m/s)
subroutine canopy_aero_ddep_katul2010(nlev, z, h_can, lad, u, d_p, rho_p, T, P, vdep_aero)
    integer, intent(in) :: nlev
    real(kind=8), intent(in) :: z(nlev), h_can, lad(nlev), u(nlev), d_p, rho_p, T(nlev), P(nlev)
    real(kind=8), intent(out) :: vdep_aero(nlev)

    ! Physical constants
    real(kind=8), parameter :: kB = 1.380649e-23   ! Boltzmann constant (J/K)
    real(kind=8), parameter :: mu_air = 1.8e-5     ! Dynamic viscosity of air (kg/m/s)
    real(kind=8), parameter :: g = 9.81            ! Gravity (m/s^2)

    integer :: i
    real(kind=8) :: D_p, Cc, V_s, Re_p, Sc, St, r_lam, r_imp, r_int, r_total
    real(kind=8) :: D_air, rho_air, nu_air

    do i = 1, nlev
        ! Calculate air density
        rho_air = P(i) / (287.05 * T(i))
        nu_air = mu_air / rho_air

        ! Cunningham slip correction factor
        Cc = 1.0 + 2.52e-7/d_p

        ! Brownian diffusion coefficient (m^2/s)
        D_air = kB * T(i) * Cc / (3.0 * pi * mu_air * d_p)

        ! Particle settling velocity (m/s)
        V_s = rho_p * g * d_p**2 * Cc / (18.0 * mu_air)

        ! Particle Reynolds number
        Re_p = u(i) * d_p / nu_air

        ! Schmidt number
        Sc = nu_air / D_air

        ! Stokes number
        St = V_s / u(i)

        ! Laminar (Brownian) resistance
        r_lam = 1.0 / (0.01 + 0.74 * D_air**0.67 * lad(i))

        ! Impaction resistance
        r_imp = 1.0 / (0.24 * St**0.6 * lad(i))

        ! Interception resistance
        r_int = 1.0 / (0.6 * d_p * lad(i))

        ! Total resistance (simplified sum)
        r_total = r_lam + r_imp + r_int

        ! Deposition velocity (m/s)
        vdep_aero(i) = 1.0 / r_total
    end do

end subroutine canopy_aero_ddep_katul2010

end module canopy_aero_ddep_mod
