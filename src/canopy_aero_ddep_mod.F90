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
    public :: canopy_aero_ddep_pleim2022

contains

!> \brief Aerosol dry deposition velocity for non-vegetated/urban areas using Pleim et al. (2022)
!> \details Implements Pleim et al. (2022) equations for urban and non-vegetated surfaces.
!>          https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022MS003050
!> \param u Array of wind speed (m/s)
!> \param d_p Aerosol particle diameter (m)
!> \param rho_p Particle density (kg/m^3)
!> \param T Air temperature (K)
!> \param P Air pressure (Pa)
!> \param surface_type Integer code for surface type (1=urban, 2=bare soil, etc.)
!> \param vdep Output: aerosol deposition velocity for surface (m/s)
    subroutine canopy_aero_ddep_pleim2022(u, d_p, rho_p, T, P, surface_type, vdep)

        use canopy_const_mod                          !< Constants for canopy models

        integer, intent(in) :: surface_type
        real(rk), intent(in) :: u, d_p, rho_p, T, P
        real(rk), intent(out) :: vdep

        ! Physical constants
        real(rk), parameter :: kB = 1.380649e-23_rk   ! Boltzmann constant (J/K)
        real(rk), parameter :: mu_air = 1.8e-5_rk     ! Dynamic viscosity of air (kg/m/s)
        real(rk), parameter :: g = 9.81_rk            ! Gravity (m/s^2)


        real(rk) :: D_air, rho_air, nu_air, Cc, V_s, Re_p, Sc, St
        real(rk) :: r_aero, v_min, v_max

        ! Calculate air density
        rho_air = P / (287.05_rk * T)
        nu_air = mu_air / rho_air

        ! Cunningham slip correction factor
        Cc = 1.0_rk + 2.52e-7_rk/d_p

        ! Brownian diffusion coefficient (m^2/s)
        D_air = kB * T * Cc / (3.0_rk * pi * mu_air * d_p)

        ! Particle settling velocity (m/s)
        V_s = rho_p * g * d_p**2.0_rk * Cc / (18.0_rk * mu_air)

        ! Particle Reynolds number
        Re_p = u * d_p / nu_air

        ! Schmidt number
        Sc = nu_air / D_air

        ! Stokes number
        St = V_s / u

        ! Pleim et al. (2022) urban/bare soil resistance parameterization
        select case (surface_type)
          case (1) ! Urban
            v_min = 0.001_rk
            v_max = 0.02_rk
            r_aero = 1.0_rk / (v_min + (v_max-v_min)*exp(-St*10.0_rk))
          case (2) ! Bare soil
            v_min = 0.0005_rk
            v_max = 0.015_rk
            r_aero = 1.0_rk / (v_min + (v_max-v_min)*exp(-St*8.0_rk))
          case default
            v_min = 0.0001_rk
            v_max = 0.01_rk
            r_aero = 1.0_rk / (v_min + (v_max-v_min)*exp(-St*5.0_rk))
        end select

        ! Total deposition velocity (m/s)
        vdep = 1.0_rk / r_aero + V_s
        ! Convert to (cm/s)
        vdep = vdep*100.0_rk

    end subroutine canopy_aero_ddep_pleim2022

!> \brief Sub-canopy aerosol dry deposition velocity following Katul et al. (2010)
!> \param nlev Number of canopy layers
!> \param lad Leaf area density profile (m^2/m^3)
!> \param z canopy model level heights (m)
!> \param hc canopy top level heights (m)
!> \param u Array of wind speed profile (m/s)
!> \param d_p Aerosol particle diameter (m)
!> \param rho_p Particle density (kg/m^3)
!> \param T Air temperature profile (K)
!> \param P Air pressure profile (Pa)
!> \param vdep Output: aerosol deposition velocity profile (m/s)
    subroutine canopy_aero_ddep_katul2010(nlev, z, hc, lad, u, d_p, rho_p, T, P, vdep)

        use canopy_const_mod                !< Constants for canopy models

        integer, intent(in) :: nlev
        real(rk), intent(in) :: lad(:), z(:), u(:), hc, d_p, rho_p, T(:), P(:)
        real(rk), intent(out) :: vdep(:)

        ! Physical constants
        real(rk), parameter :: kB = 1.380649e-23_rk   ! Boltzmann constant (J/K)
        real(rk), parameter :: mu_air = 1.8e-5_rk     ! Dynamic viscosity of air (kg/m/s)
        real(rk), parameter :: g = 9.81_rk            ! Gravity (m/s^2)

        integer :: i
        real(rk) :: Cc, V_s, Re_p, Sc, St, r_lam, r_imp, r_int, r_total
        real(rk) :: D_air, rho_air, nu_air

        do i = 1, nlev
            if (z(i) .gt. 0.0 .and. z(i) .le. hc) then  !< Above ground level and at/below canopy top
                ! Calculate air density
                rho_air = (P(i)*100.0_rk) / (287.05_rk * T(i))  !need to convert pressure profile from mb to Pa
                nu_air = mu_air / rho_air

                ! Cunningham slip correction factor
                Cc = 1.0_rk + 2.52e-7_rk/d_p

                ! Brownian diffusion coefficient (m^2/s)
                D_air = kB * T(i) * Cc / (3.0_rk * pi * mu_air * d_p)

                ! Particle settling velocity (m/s)
                V_s = rho_p * g * d_p**2.0_rk * Cc / (18.0_rk * mu_air)

                ! Particle Reynolds number
                Re_p = u(i) * d_p / nu_air

                ! Schmidt number
                Sc = nu_air / D_air

                ! Stokes number
                St = V_s / u(i)

                ! Laminar (Brownian) resistance
                r_lam = 1.0_rk / (0.01_rk + 0.74_rk * D_air**0.67_rk * lad(i))

                ! Impaction resistance
                r_imp = 1.0_rk / (0.24_rk * St**0.6_rk * lad(i))

                ! Interception resistance
                r_int = 1.0_rk / (0.6_rk * d_p * lad(i))

                ! Total resistance (simplified sum)
                r_total = r_lam + r_imp + r_int

                ! Deposition velocity (m/s)
                vdep(i) = 1.0_rk / r_total

                ! Convert to (cm/s)
                vdep(i) = vdep(i)*100.0_rk

            else
                vdep(i) = 0.0_rk
            end if
        end do

    end subroutine canopy_aero_ddep_katul2010

end module canopy_aero_ddep_mod
