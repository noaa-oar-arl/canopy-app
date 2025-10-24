

# File canopy\_aero\_ddep\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_aero\_ddep\_mod.F90**](canopy__aero__ddep__mod_8F90.md)

[Go to the documentation of this file](canopy__aero__ddep__mod_8F90.md)


```Fortran


module canopy_aero_ddep_mod

    implicit none
    private
    public :: canopy_aero_ddep_subveg
    public :: canopy_aero_ddep_pleim2022

contains

    subroutine canopy_aero_ddep_pleim2022(ustar, Ra, u, d_p, rho_p, T, P, surface_type, vdep)

        use canopy_const_mod                          !< constants for canopy models

        integer, intent(in) :: surface_type
        real(rk), intent(in) :: ustar, ra, u, d_p, rho_p, t, p
        real(rk), intent(out) :: vdep

        ! Physical constants
        real(rk), parameter :: kb = 1.380649e-23_rk   ! Boltzmann constant (J/K)
        real(rk), parameter :: mu_air = 1.8e-5_rk     ! Dynamic viscosity of air (kg/m/s)
        real(rk), parameter :: g = 9.81_rk            ! Gravity (m/s^2)


        real(rk) :: d_air, rho_air, nu_air, cc, v_s, re_p, sc, st
        real(rk) :: r_aero, eb, eim, aa, bb, fwc

        ! Calculate air density
        rho_air = p / (287.05_rk * t)
        nu_air = mu_air / rho_air

        ! Cunningham slip correction factor
        cc = 1.0_rk + 2.52e-7_rk/d_p

        ! Brownian diffusion coefficient (m^2/s)
        d_air = kb * t * cc / (3.0_rk * pi * mu_air * d_p)

        ! Particle settling velocity (m/s)
        v_s = rho_p * g * d_p**2.0_rk * cc / (18.0_rk * mu_air)

        ! Particle Reynolds number
        re_p = u * d_p / nu_air

        ! Schmidt number
        sc = nu_air / d_air

        ! Stokes number
        st = v_s / u

        ! Pleim et al. (2022) urban/bare soil resistance parameterizations for smooth surfaces
        ! Big-leaf approach outside vegetative canopies
        select case (surface_type)
          case (1) ! Urban
            ! Laminar (Brownian) term
            eb=(1.0_rk/3.0_rk) * sc**(-2.0_rk/3.0_rk)
            ! Impaction term
            eim = 10.**(-3.0_rk/st)
            !Urban resistance
            r_aero = 1.0_rk /(2.0_rk*ustar * (eb + eim)) !Assume BAI = 2 for developed areas
          case (2) ! Bare soil
            ! Laminar (Brownian) term
            eb=(1.0_rk/3.0_rk) * sc**(-2.0_rk/3.0_rk)
            ! Impaction term
            eim = 10.0_rk**(-3.0_rk/st)
            !Bare ground/soil resistance
            r_aero = 1.0_rk / (ustar * (eb + eim))
          case default !Other Water Surfaces
            !additional terms due to wave breaking
            aa = 8.46e-5_rk + (1.63e-6_rk*(t-273.15)) + (-3.35e-8_rk*(t-273.15)**2.0)
            bb = 3.354 + (-0.062*(t-273.15))
            fwc = aa*(bb+u)**2.0
            ! Laminar (Brownian) term
            eb=(1.0_rk - fwc)*(1.0_rk/3.0_rk) * sc**(-2.0_rk/3.0_rk) + fwc*(ustar/u)
            ! Impaction term
            eim = 10.0_rk**(-3.0_rk/st)
            !Water resistance
            r_aero = 1.0_rk / (ustar * (eb + eim))
        end select

        ! Total deposition velocity (m/s)
        vdep = v_s / (1.0_rk - exp(-1.0_rk*v_s*(r_aero+ra)))
        ! Convert to (cm/s)
        vdep = vdep*100.0_rk

    end subroutine canopy_aero_ddep_pleim2022

    subroutine canopy_aero_ddep_subveg(nlev, z, hc, lad, u, d_p, rho_p, T, P, vdep_opt, Ra, modres, ustar, &
        vtype, vdep)

        use canopy_const_mod                !< constants for canopy models

        integer, intent(in) :: nlev, vdep_opt, vtype
        real(rk), intent(in) :: lad(:), z(:), u(:), hc, d_p, rho_p, t(:), p(:), ra, modres, ustar
        real(rk), intent(out) :: vdep(:)

        ! Physical constants
        real(rk), parameter :: kb = 1.380649e-23_rk   ! Boltzmann constant (J/K)
        real(rk), parameter :: mu_air = 1.8e-5_rk     ! Dynamic viscosity of air (kg/m/s)
        real(rk), parameter :: g = 9.81_rk            ! Gravity (m/s^2)

        integer :: i
        real(rk) :: cc, v_s, re_p, sc, stl, sth, eb, eim, r_int, r_aero, r_total
        real(rk) :: d_air, rho_air, nu_air, laix, al, ah, fmicro

        do i = 1, nlev
            if (z(i) .gt. 0.0 .and. z(i) .le. hc) then  
                ! Calculate air density
                rho_air = (p(i)*100.0_rk) / (287.05_rk * t(i))  !need to convert pressure profile from mb to Pa
                nu_air = mu_air / rho_air

                ! Cunningham slip correction factor
                cc = 1.0_rk + 2.52e-7_rk/d_p

                ! Brownian diffusion coefficient (m^2/s)
                d_air = kb * t(i) * cc / (3.0_rk * pi * mu_air * d_p)

                ! Particle settling velocity (m/s)
                v_s = rho_p * g * d_p**2.0_rk * cc / (18.0_rk * mu_air)

                ! Particle Reynolds number
                re_p = u(i) * d_p / nu_air

                ! Schmidt number
                sc = nu_air / d_air

                ! Stokes number
                !St = V_s / u(i)


                !Calculate resistances based on Pleim et al. (2022)
                !Including microscale leaf effects
                !Note:  Particle rebound effects are also not included (R = 1)
                if (vtype .ge. 1 .and. vtype .le. 2) then !needleleaf
                    al = 2.0_rk/1000.0_rk  !2 mm    (converted to meters)
                    ah = 0.5_rk/1.0e6_rk   !0.5 um  (converted to meters)
                    fmicro = 0.008_rk
                elseif (vtype .ge. 3 .and. vtype .le. 5) then !broadleaf
                    al = 10.0_rk/1000.0_rk !10 mm   (converted to meters)
                    ah = 1.0_rk/1.0e6_rk   !0.5 um  (converted to meters)
                    fmicro = 0.008_rk
                else !other vegetated canopies, e.g., grasslands, shrublands, etc.
                    al = 0.5_rk/1000.0_rk  !0.5 mm  (converted to meters)
                    ah = 0.5_rk/1.0e6_rk   !0.5 um  (converted to meters)
                    fmicro = 0.002_rk
                end if

                ! Stokes number for vegetated surfaces (Pleim et al., 2022)
                stl = v_s*ustar/g*al  !macro/leaf scale
                sth = v_s*ustar/g*ah  !micro/hair sclae

                ! Laminar (Brownian) term
                eb=(1.0_rk/3.0_rk) * sc**(-2.0_rk/3.0_rk)
                ! Impaction term
                !New Impaction term following Pleim et al. (2022)
                eim = (1.0_rk - fmicro)*(stl**2.0_rk)/(1.0_rk + stl**2.0_rk) + &
                    fmicro*(sth**2.0_rk)/(1.0_rk + sth**2.0_rk)
                !Canopy layer fractional LAI
                laix=lad(i)*modres
                !Resistance
                r_aero = 1.0_rk / (laix * ustar* (eb + eim))

                !Add interception resistance based on Katul, Petroff, or Zhang
                if (vdep_opt == 0) then
                    ! Katul Interception resistance scaled to LAD
                    r_int = 1.0_rk / (0.6_rk * d_p * lad(i))
                else if (vdep_opt == 1) then
                    ! Petroff Interception resistance scaled to LAD
                    r_int = 1.0_rk / (0.5_rk * d_p * lad(i))
                else
                    ! Zhang Interception resistance scaled to LAD
                    r_int = 1.0_rk / (0.5_rk * d_p * lad(i))
                end if

                ! Total aerosol resistance (parallel combination)
                r_total = 1.0_rk / (1.0_rk/r_aero + 1.0_rk/r_int)

                ! Total deposition velocity (m/s)
                vdep(i) = v_s / (1.0_rk - exp(-1.0_rk*v_s*(r_total+ra)))
                ! Convert to (cm/s)
                vdep(i) = vdep(i)*100.0_rk

            else
                vdep(i) = 0.0_rk
            end if
        end do

    end subroutine canopy_aero_ddep_subveg

end module canopy_aero_ddep_mod
```


