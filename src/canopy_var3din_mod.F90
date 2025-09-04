!> \file canopy_var3din_mod.F90
!! \brief 3D Variable Input Module
!! \details This module contains subroutines for processing 3D variable inputs,
!! particularly for converting GEDI PAVD (Plant Area Volume Density) profiles
!! into fractional foliage shape functions. The module handles interpolation
!! of observed PAVD profiles to user-defined canopy model resolutions.
!!
!! \author Patrick C. Campbell
!! \date July 2023
!!
!! \references
!! Massman, W.J., Forthofer, J.M., and Finney, M.A.: An improved
!! canopy wind model for predicting wind adjustment factors
!! and wildland fire behavior. Canadian Journal of Forest Research.
!! 47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354

!> \defgroup var3din_mod 3D Variable Input Module
!! \brief Module for processing 3D canopy structure inputs
!! \{

module canopy_var3din_mod

    implicit none

contains

!> \brief Compute integral of incremental fractional foliage shape function from GEDI PAVD
!! \details This subroutine converts interpolated GEDI 3D PAVD (Plant Area Volume Density)
!! profiles into fractional foliage shape functions using the algorithms from Massman et al. (2017).
!! The process includes:
!! - Interpolating input PAVD data to the canopy model vertical resolution
!! - Determining the height of maximum foliage area density (ZCANMAX) from observed PAVD
!! - Calculating incremental foliage shape functions using Gaussian distributions
!! - Computing fractional cumulative foliage distributions
!! - Integrating the foliage shape functions for canopy structure parameterization
!!
!! \param[in] ZCANMAX_IN Input height of maximum foliage area density (z/h) (nondimensional)
!! \param[in] SIGMAU Standard deviation of shape function above zcanmax (z/h)
!! \param[in] SIGMA1 Standard deviation of shape function below zcanmax (z/h)
!! \param[in] FCH Grid cell canopy height (m) from GEDI
!! \param[in] ZHC Dimensionless height coordinate (z/h)
!! \param[in] PAVD_IN Plant Area Volume Density profile (m²/m³)
!! \param[in] PAVD_LEVS Associated mid-level heights for PAVD data (m)
!! \param[out] FAFRACZINT Integral of incremental fractional foliage shape function
    SUBROUTINE CANOPY_PAVD2FAFRAC ( ZCANMAX_IN, SIGMAU, SIGMA1, FCH, &
        ZHC, PAVD_IN, PAVD_LEVS, FAFRACZINT )
        use canopy_const_mod, ONLY: rk                                           !> constants for canopy models
        use canopy_utils_mod, ONLY: interp_linear1_internal,IntegrateTrapezoid   !> utilities for canopy models

!> \defgroup var3din_inputs Input Variables
!! \brief Input parameters for PAVD processing
!! \{
        REAL(RK),    INTENT( IN )  :: FCH                   !> Grid cell canopy height (m) from GEDI
        REAL(RK),    INTENT( IN )  :: ZCANMAX_IN            !> Input height of maximum foliage area density (z/h) (nondimensional)
        REAL(RK),    INTENT( IN )  :: SIGMAU                !> Standard deviation of shape function above zcanmax (z/h)
        REAL(RK),    INTENT( IN )  :: SIGMA1                !> Standard deviation of shape function below zcanmax (z/h)
        REAL(RK),    INTENT( IN )  :: ZHC(:)                !> z/h (dimensionless)
        REAL(RK),    INTENT( IN )  :: PAVD_IN(:)            !> Plant Area Volume Density (PAVD) profile (m2/m3)
        REAL(RK),    INTENT( IN )  :: PAVD_LEVS(:)          !> Associated mid-level heights associated with PAVD (m)
!> \}

!> \defgroup var3din_outputs Output Variables
!! \brief Output parameters for foliage shape calculations
!! \{
        REAL(RK),    INTENT( OUT ) :: FAFRACZINT(:)         !> integral of incremental fractional foliage shape function
!> \}

!> \defgroup var3din_local_vars Local Variables
!! \brief Local variables for PAVD processing
!! \{
        INTEGER                     :: i, lev               !> Loop counters
        REAL(RK)                    :: ZK(SIZE(ZHC))        !> Actual height coordinates (m)
        REAL(RK)                    :: PAVD_INTERP(SIZE(ZHC)) !> Interpolated PAVD profile (m²/m³)

        REAL(RK)                    :: ZCANMAX               !> Height of maximum foliage area density (z/h) (nondimensional)
        REAL(RK), allocatable       :: fainc(:)              !> incremental foliage shape function
        REAL(RK), allocatable       :: fafracz(:)            !> incremental fractional foliage shape function
        REAL(RK)                    :: fatot                 !> integral of total fractional foliage shape function
!> \}

        !> \brief Convert dimensionless heights to actual heights
        ZK = ZHC*FCH

        !> \brief Interpolate input PAVD at set levels to user desired canopy model resolution
        PAVD_INTERP = 0.0_rk  !Initialize PAVD_INTERP = 0
        do lev=1, SIZE(PAVD_LEVS) - 1
            do i=2, SIZE(ZK)  !loop over only levels ABOVE ground
                if (ZK(i) .le.  PAVD_LEVS(1)) then
                    PAVD_INTERP(i)   = PAVD_IN(1)
                end if
                if (ZK(i) .ge.  PAVD_LEVS(lev) .and. ZK(i) .le.  PAVD_LEVS(lev+1)) then
                    PAVD_INTERP(i)   = interp_linear1_internal((/ PAVD_LEVS(lev),PAVD_LEVS(lev+1) /), &
                        (/ PAVD_IN(lev),PAVD_IN(lev+1) /),ZK(i))
                end if
            end do
        end do

        !> \brief Initialize canopy profile dependent variables
        ZCANMAX             = 0.0_rk

        !> \brief Bottom-up loop to determine ZCANMAX from observed PAVD
        !! \details Find the height where PAVD reaches its maximum value
        do i=2, SIZE(ZK)
            if (PAVD_INTERP(i) .ge. maxval(PAVD_INTERP) ) then
                ZCANMAX = ZK(i)/FCH
                if(ZCANMAX .gt. 1.0_rk) then !if ZK (at max GEDI PAVD height) > GEDI FCH (inconsistent!)
                    ZCANMAX = ZCANMAX_IN       !override with Massman Input ZCANMAX
                end if
                exit
            end if
        end do

        !> \brief Calculate foliage shape functions using Massman et al. (2017) Eqs. 1-2
        if(.not.allocated(fainc))      allocate(fainc(SIZE(ZK)))
        if(.not.allocated(fafracz))    allocate(fafracz(SIZE(ZK)))

        !> \brief Initialize canopy profile dependent variables
        fainc             = 0.0_rk
        fafracz           = 0.0_rk

        !> \brief Calculate canopy/foliage distribution shape profile using Gaussian functions
        !! \details Uses different standard deviations above and below ZCANMAX
        do i=1, SIZE(ZK)
            if (ZHC(i) >= ZCANMAX .and. ZHC(i) <= 1.0) then
                fainc(i) = exp((-1.0*((ZHC(i)-ZCANMAX)**2.0))/SIGMAU**2.0)
            else if (ZHC(i) >= 0.0 .and. ZHC(i) <= ZCANMAX) then
                fainc(i) = exp((-1.0*((ZCANMAX-ZHC(i))**2.0))/SIGMA1**2.0)
            end if
        end do
        fatot = IntegrateTrapezoid(ZHC,fainc)

        !> \brief Calculate normalized plant distribution function and its integral
        !! \details Normalize by total integral and compute cumulative distribution
        do i=1, SIZE(ZK)
            fafracz(i) = fainc(i)/fatot
            FAFRACZINT(i) = IntegrateTrapezoid(ZHC(1:i),fafracz(1:i))
        end do

    END SUBROUTINE CANOPY_PAVD2FAFRAC

!> \}

end module canopy_var3din_mod
