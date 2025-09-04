!> \file canopy_const_mod.F90
!! \brief Constants Module
!! \details This module contains fundamental physical and mathematical constants
!! for air quality and canopy modeling applications. Constants are derived from
!! standard references including the CRC Handbook and atmospheric sciences texts.
!!
!! \author Patrick C. Campbell
!! \date July 2022
!!
!! \references
!! - CRC76: "CRC Handbook of Chemistry and Physics (76th Ed)", CRC Press, 1995
!! - Hobbs, P.V.: "Basic Physical Chemistry for the Atmospheric Sciences",
!!   Cambridge Univ. Press, 206 pp, 1995
!! - Snyder, J.P.: "Map Projections-A Working Manual, U.S. Geological Survey
!!   Paper 1395 U.S.GPO, Washington, DC, 1987
!! - Stull, R. B.: "An Introduction to Boundary Layer Meteorology", Kluwer,
!!   Dordrecht, 1988

!> \defgroup const_mod Constants Module
!! \brief Module containing fundamental constants for canopy modeling
!! \{

MODULE canopy_const_mod

    IMPLICIT NONE

!> \defgroup precision_const Precision Constants
!! \brief Precision and fill value constants
!! \{

    !> \brief Selected real kind for double precision (15 decimal digits, 307 exponent range)
    INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
    !> \brief NetCDF fill value for missing real data
    REAL(rk),    PARAMETER :: fillreal  = -9.0e20

!> \}

!> \defgroup geom_const Geometric Constants
!! \brief Mathematical and geometric constants
!! \{

    !> \brief Pi (double precision: 3.14159265358979324)
    REAL(RK),       PARAMETER     :: pi = 3.14159265358979324_rk

    !> \brief Pi/180 conversion factor from degrees to radians [rad/deg]
    REAL(RK),          PARAMETER     :: pi180 = pi / 180.0_rk

!> \}

!> \defgroup geodetic_const Geodetic Constants
!! \brief Earth geometry and geodetic constants
!! \{

    !> \brief Radius of Earth [m]
    !! \details Radius of sphere having same surface area as Clarke ellipsoid of 1866
    !! WGS84 arithmetic mean radius (Source: Snyder, 1987)
    REAL(RK),          PARAMETER     :: rearth = 6371008.8_rk

    !> \brief Length of a sidereal day [sec] (Source: CRC76, pp. 14-6)
    REAL(RK),          PARAMETER     :: siday = 86164.09_rk

    !> \brief Mean gravitational acceleration [m/sec²]
    !! \details Mean of polar and equatorial values (Source: CRC76, pp. 14-6)
    REAL(RK),          PARAMETER     :: grav = 9.80622_rk

    !> \brief Latitude degrees to meters conversion factor [m/deg]
    REAL(RK)                         :: dg2m

    !> \brief Solar Constant [W/m²] (Source: CRC76, pp. 14-2)
    REAL(RK),          PARAMETER     :: solcnst = 1373.0

!> \}


!> \defgroup fundamental_const Fundamental Constants
!! \brief Fundamental physical constants (Source: CRC76, pp. 1-1 to 1-6)
!! \{

    !> \brief Avogadro's Constant [number/mol]
    REAL(RK),          PARAMETER     :: avo = 6.02214076d23

    !> \brief Universal gas constant [J/mol-K]
    REAL(RK),          PARAMETER     :: rgasuniv = 8.31446261815324_rk

    !> \brief Standard atmosphere [Pa]
    REAL(RK),          PARAMETER     :: stdatmpa = 101325.0

    !> \brief Standard temperature [K]
    REAL(RK),          PARAMETER     :: stdtemp = 273.15_rk

    !> \brief Stefan-Boltzmann constant [W/(m² K⁴)]
    REAL(RK),          PARAMETER     :: stfblz = 5.67037442d-8

    !> \brief Molar volume of ideal gas at STP [L/mol] (Non-MKS units)
    REAL(RK),          PARAMETER     :: molvol = 22.4139695_rk

!> \}


!> \defgroup atmos_const Atmospheric Constants
!! \brief Constants for atmospheric thermodynamics and chemistry
!! \{

    !> \brief Mean molecular weight for dry air [g/mol]
    !! \details 78.06% N2, 21% O2, and 0.943% A on a mole fraction basis
    !! (Source: Hobbs, 1995, pp. 69-70)
    REAL(RK),          PARAMETER     :: mwair = 28.9628_rk

    !> \brief Dry-air gas constant [J/kg-K]
    !! \details 287.07548994 J/kg-K
    REAL(RK),          PARAMETER     :: rdgas = 1.0e3 * rgasuniv / mwair

    !> \brief Mean molecular weight for water vapor [g/mol]
    REAL(RK),          PARAMETER     :: mwwat = 18.01528_rk

    !> \brief Gas constant for water vapor [J/kg-K]
    !! \details 461.52492604 J/kg-K
    REAL(RK),          PARAMETER     :: rwvap = 1.0e3 * rgasuniv / mwwat

    !> \brief Specific heat of dry air at constant pressure [J/kg-K]
    !! \details 1004.7642148 J/kg-K. Calculated assuming dry air is classical ideal gas
    REAL(RK),          PARAMETER     :: cpd = 7.0 * rdgas / 2.0

    !> \brief Specific heat of dry air at constant volume [J/kg-K]
    !! \details 717.68872485 J/kg-K. Calculated assuming dry air is classical ideal gas
    REAL(RK),          PARAMETER     :: cvd = 5.0 * rdgas / 2.0

    !> \brief Specific heat for water vapor at constant pressure [J/kg-K]
    !! \details 1846.0997042 J/kg-K. Calculated assuming water vapor is classical ideal gas
    REAL(RK),          PARAMETER     :: cpwvap = 4.0 * rwvap

    !> \brief Specific heat for water vapor at constant volume [J/kg-K]
    !! \details 1384.5747781 J/kg-K. Calculated assuming water vapor is classical ideal gas
    REAL(RK),          PARAMETER     :: cvwvap = 3.0 * rwvap

    !> \brief Vapor pressure of water at 0°C [Pa] (Source: CRC76 pp. 6-15)
    REAL(RK),          PARAMETER     :: vp0 = 611.29_rk

    !> \brief Latent heat of vaporization of water at 0°C [J/kg]
    !! \details Values from p. 641 of Stull (1988)
    REAL(RK),          PARAMETER     :: lv0 = 2.501e6

    !> \brief Rate of change of latent heat of vaporization w.r.t. temperature [J/kg-K]
    REAL(RK),          PARAMETER     :: dlvdt = 2370.0

    !> \brief Latent heat of fusion of water at 0°C [J/kg]
    REAL(RK),          PARAMETER     :: lf0 = 3.34e5

    !> \brief Von Karman constant (dimensionless)
    REAL(RK),          PARAMETER     :: vonk = 0.4_rk

    !> \brief Stability parameter beta for neutral conditions (dimensionless)
    !! \details Source: Bonan et al. (2018) https://doi.org/10.5194/gmd-11-1467-2018
    REAL(RK),          PARAMETER     :: beta_n = 0.35_rk

!> \}

!> \defgroup bioemi_const Biogenic Emission Constants
!! \brief Constants for biogenic emission calculations
!! \{

    !> \brief E-folding time for long-term past conditions for biogenic emissions [days]
    REAL(RK), PARAMETER  :: tau_days  = 5.0_rk

    !> \brief E-folding time for short-term past conditions for biogenic emissions [hours]
    REAL(RK), PARAMETER  :: tau_hours = 12.0_rk

!> \}

!> \}

END MODULE canopy_const_mod
