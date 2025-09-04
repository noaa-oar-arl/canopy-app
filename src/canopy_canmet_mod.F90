!> \file canopy_canmet_mod.F90
!> \brief Canopy meteorological and surface input variable definitions
!> \details This module contains type definitions and variables for meteorological
!>          and surface input data used by the canopy model, including 2D and 3D
!>          input variables, canopy profile data, and reference conditions.
!> \author P. C. Campbell
!> \date 03 Oct 2022

!> \defgroup canopy_met_input Canopy Meteorological Input Variables
!> \brief Meteorological and surface input variables for canopy model

MODULE canopy_canmet_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Met/Sfc Input Variable Descriptions
! Purpose:  Contains canopy met and sfc input  variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk

    IMPLICIT NONE
!! .... defines canopy options (read from user namelist)

    !> \ingroup canopy_met_input
    !> \{

    !> \brief 2D meteorological and surface input variables
    !> \details Generic 2D met/sfc input variables that should be passed to canopy calculations
    TYPE :: variable_type
        !> \brief Latitude coordinate
        !> \details Latitude of cell/point
        !! \param units degrees
        real(rk)   :: lat

        !> \brief Longitude coordinate
        !> \details Longitude of cell/point
        !! \param units degrees
        real(rk)   :: lon

        !> \brief Canopy height
        !> \details Canopy height
        !! \param units meters (m)
        real(rk)   :: ch

        !> \brief U-component wind speed
        !> \details U wind speed at reference height above canopy
        !! \param units m/s
        real(rk)   :: ugrd10m

        !> \brief V-component wind speed
        !> \details V wind speed at reference height above canopy
        !! \param units m/s
        real(rk)   :: vgrd10m

        !> \brief Clumping index
        !> \details Clumping index for canopy structure
        !! \param units dimensionless
        real(rk)   :: clu

        !> \brief Leaf area index
        !> \details Leaf area index
        !! \param units m²/m²
        real(rk)   :: lai

        !> \brief Vegetation type
        !> \details Vegetation type classification
        integer    :: vtype

        !> \brief Canopy fraction
        !> \details Canopy fraction of grid cell
        !! \param units dimensionless fraction
        real(rk)   :: canfrac

        !> \brief Friction velocity
        !> \details Friction velocity (u*)
        !! \param units m/s
        real(rk)   :: fricv

        !> \brief Cosine of solar zenith angle
        !> \details Cosine of solar zenith angle
        !! \param units dimensionless
        real(rk)   :: csz

        !> \brief Surface roughness length
        !> \details Surface roughness length
        !! \param units meters (m)
        real(rk)   :: sfcr

        !> \brief Monin-Obukhov length
        !> \details Monin-Obukhov length
        !! \param units meters (m)
        real(rk)   :: mol

        !> \brief Fire radiative power
        !> \details Fire radiative power
        !! \param units MW
        real(rk)   :: frp

        !> \brief Reference height above canopy
        !> \details Reference height above the canopy
        !! \param units meters (m)
        real(rk)   :: href

        !> \brief Soil type
        !> \details Soil type classification
        integer    :: sotyp

        !> \brief Surface pressure
        !> \details Surface pressure
        !! \param units hPa
        real(rk)   :: pressfc

        !> \brief Downward shortwave radiation
        !> \details Instantaneous downward shortwave radiation
        !! \param units W/m²
        real(rk)   :: dswrf

        !> \brief Surface sensible heat flux
        !> \details Instantaneous surface sensible heat net flux
        !! \param units W/m²
        real(rk)   :: shtfl

        !> \brief Surface temperature
        !> \details Surface temperature
        !! \param units K
        real(rk)   :: tmpsfc

        !> \brief 2-meter temperature
        !> \details 2-meter temperature
        !! \param units K
        real(rk)   :: tmp2m

        !> \brief 2-meter specific humidity
        !> \details 2-meter specific humidity
        !! \param units kg/kg
        real(rk)   :: spfh2m

        !> \brief Planetary boundary layer height
        !> \details Height of planetary boundary layer
        !! \param units meters (m)
        real(rk)   :: hpbl

        !> \brief Precipitation rate
        !> \details Mass precipitation rate
        !! \param units kg/m²/s
        real(rk)   :: prate_ave

        !> \brief Soil moisture level 1
        !> \details Volumetric soil moisture level 1
        !! \param units m³/m³
        real(rk)   :: soilw1

        !> \brief Soil moisture level 2
        !> \details Volumetric soil moisture level 2
        !! \param units m³/m³
        real(rk)   :: soilw2

        !> \brief Soil moisture level 3
        !> \details Volumetric soil moisture level 3
        !! \param units m³/m³
        real(rk)   :: soilw3

        !> \brief Soil moisture level 4
        !> \details Volumetric soil moisture level 4
        !! \param units m³/m³
        real(rk)   :: soilw4

        !> \brief Wilting point
        !> \details Wilting point
        !! \param units proportion
        real(rk)   :: wilt

        !> \brief Ozone W126 values
        !> \details Ozone W126 values
        !! \param units ppm-hours
        real(rk)   :: ozone_w126

        !> \brief Soil temperature level 1
        !> \details Soil temperature level 1
        !! \param units K
        real(rk)   :: soilt1

        !> \brief Soil temperature level 2
        !> \details Soil temperature level 2
        !! \param units K
        real(rk)   :: soilt2

        !> \brief Soil temperature level 3
        !> \details Soil temperature level 3
        !! \param units K
        real(rk)   :: soilt3

        !> \brief Soil temperature level 4
        !> \details Soil temperature level 4
        !! \param units K
        real(rk)   :: soilt4

        !> \brief First model layer air temperature
        !> \details 1st model layer air temperature above ground
        !! \param units K
        real(rk)   :: tmp_hyblev1

        !> \brief Average ground snow cover
        !> \details Average percent ground snow cover
        !! \param units percent (%)
        real(rk)   :: snowc_ave
        !> \brief Average ground or water ice cover
        !> \details Average fraction ground or water ice cover
        !! \param units dimensionless fraction
        real(rk)   :: icec
    end TYPE variable_type

    !> \brief Allocated array for 1D meteorological variables
    !> \details Array for storing 1D meteorological input variables
    type(variable_type), allocatable :: variables( : ), variables_2d( : , :)

    !> \brief 1D vertical level input variables
    !> \details Generic 3D input variables for vertical levels
    TYPE :: variable_type_1d
        !> \brief Input mid-level heights
        !> \details Input mid-level heights associated with 3D input option
        !! \param units meters (m)
        real(rk)   :: lev
    end TYPE variable_type_1d

    !> \brief Allocated array for 1D vertical variables
    !> \details Array for storing 1D vertical level variables
    type(variable_type_1d), allocatable :: variables_1d( : )

    !> \brief 3D Plant Area Volume Density variables
    !> \details Type for 3D PAVD input variables
    TYPE :: variable_type_3d
        !> \brief Plant Area Volume Density profile
        !> \details Plant Area Volume Density (PAVD) profile
        !! \param units m²/m³
        real(rk)   :: pavd
    end TYPE variable_type_3d

    !> \brief Allocated array for 3D PAVD variables
    !> \details Array for storing 3D PAVD variables
    type(variable_type_3d), allocatable :: variables_3d( : , : , :)

    !> \brief Canopy profile input variables with 14 levels
    !> \details Generic set of observed canopy profile input variable levels (14) from point text file
    TYPE :: variable_type_can
        !> \brief Latitude coordinate
        !> \details Latitude of cell/point
        !! \param units degrees
        real(rk)   :: lat

        !> \brief Longitude coordinate
        !> \details Longitude of cell/point
        !! \param units degrees
        real(rk)   :: lon

        !> \brief Input canopy profile level 1
        !> \details Input canopy profile level 1
        !! \param units meters (m)
        real(rk)   :: lev01

        !> \brief Input canopy PAVD profile level 1
        !> \details Input canopy PAVD profile level 1
        !! \param units m²/m³
        real(rk)   :: pavd01

        !> \brief Input canopy profile level 2
        !! \param units meters (m)
        real(rk)   :: lev02
        !> \brief Input canopy PAVD profile level 2
        !! \param units m²/m³
        real(rk)   :: pavd02
        !> \brief Input canopy profile level 3
        !! \param units meters (m)
        real(rk)   :: lev03
        !> \brief Input canopy PAVD profile level 3
        !! \param units m²/m³
        real(rk)   :: pavd03
        !> \brief Input canopy profile level 4
        !! \param units meters (m)
        real(rk)   :: lev04
        !> \brief Input canopy PAVD profile level 4
        !! \param units m²/m³
        real(rk)   :: pavd04
        !> \brief Input canopy profile level 5
        !! \param units meters (m)
        real(rk)   :: lev05
        !> \brief Input canopy PAVD profile level 5
        !! \param units m²/m³
        real(rk)   :: pavd05
        !> \brief Input canopy profile level 6
        !! \param units meters (m)
        real(rk)   :: lev06
        !> \brief Input canopy PAVD profile level 6
        !! \param units m²/m³
        real(rk)   :: pavd06
        !> \brief Input canopy profile level 7
        !! \param units meters (m)
        real(rk)   :: lev07
        !> \brief Input canopy PAVD profile level 7
        !! \param units m²/m³
        real(rk)   :: pavd07
        !> \brief Input canopy profile level 8
        !! \param units meters (m)
        real(rk)   :: lev08
        !> \brief Input canopy PAVD profile level 8
        !! \param units m²/m³
        real(rk)   :: pavd08
        !> \brief Input canopy profile level 9
        !! \param units meters (m)
        real(rk)   :: lev09
        !> \brief Input canopy PAVD profile level 9
        !! \param units m²/m³
        real(rk)   :: pavd09
        !> \brief Input canopy profile level 10
        !! \param units meters (m)
        real(rk)   :: lev10
        !> \brief Input canopy PAVD profile level 10
        !! \param units m²/m³
        real(rk)   :: pavd10
        !> \brief Input canopy profile level 11
        !! \param units meters (m)
        real(rk)   :: lev11
        !> \brief Input canopy PAVD profile level 11
        !! \param units m²/m³
        real(rk)   :: pavd11
        !> \brief Input canopy profile level 12
        !! \param units meters (m)
        real(rk)   :: lev12
        !> \brief Input canopy PAVD profile level 12
        !! \param units m²/m³
        real(rk)   :: pavd12
        !> \brief Input canopy profile level 13
        !! \param units meters (m)
        real(rk)   :: lev13
        !> \brief Input canopy PAVD profile level 13
        !! \param units m²/m³
        real(rk)   :: pavd13
        !> \brief Input canopy profile level 14
        !! \param units meters (m)
        real(rk)   :: lev14
        !> \brief Input canopy PAVD profile level 14
        !! \param units m²/m³
        real(rk)   :: pavd14
    end TYPE variable_type_can

    !> \brief Allocated array for canopy profile variables
    !> \details Array for storing canopy profile variables
    type(variable_type_can), allocatable :: variables_can( : )


    !> \brief Met/Sfc variable reassignment names above reference conditions from the model
    !> \details Reference meteorological and surface variables used in canopy calculations
    !> \{

    !> \brief Reference latitude
    !> \details Latitude of cell/point
    !! \param units degrees
    real(rk)       ::    latref

    !> \brief Reference longitude
    !> \details Longitude of cell/point
    !! \param units degrees
    real(rk)       ::    lonref

    !> \brief Reference canopy height
    !> \details Input canopy height
    !! \param units meters (m)
    real(rk)       ::    hcmref

    !> \brief Reference U wind speed
    !> \details Input above canopy/reference 10-m U wind speed
    !! \param units m/s
    real(rk)       ::    uref

    !> \brief Reference V wind speed
    !> \details Input above canopy/reference 10-m V wind speed
    !! \param units m/s
    real(rk)       ::    vref

    !> \brief Reference bulk wind speed
    !> \details Input above canopy/reference 10-m model wind speed
    !! \param units m/s
    real(rk)       ::    ubzref

    !> \brief Reference clumping index
    !> \details Input canopy clumping index
    !! \param units dimensionless
    real(rk)       ::    cluref

    !> \brief Reference leaf area index
    !> \details Input leaf area index
    !! \param units m²/m²
    real(rk)       ::    lairef

    !> \brief Reference vegetation type
    !> \details Input vegetation type (VIIRS)
    integer        ::    vtyperef

    !> \brief Reference canopy fraction
    !> \details Input canopy fraction of grid cell
    !! \param units dimensionless fraction
    real(rk)       ::    canfracref

    !> \brief Reference friction velocity
    !> \details Input friction velocity
    !! \param units m/s
    real(rk)       ::    ustref

    !> \brief Reference cosine of zenith angle
    !> \details Input cosine of zenith angle
    !! \param units dimensionless
    real(rk)       ::    cszref

    !> \brief Reference surface roughness length
    !> \details Input total/surface roughness length
    !! \param units meters (m)
    real(rk)       ::    z0ref

    !> \brief Reference Monin-Obukhov length
    !> \details Input Monin-Obukhov Length
    !! \param units meters (m)
    real(rk)       ::    molref

    !> \brief Reference fire radiative power
    !> \details Input fire radiative power
    !! \param units MW
    real(rk)       ::    frpref

    !> \brief Reference height above canopy
    !> \details Reference height above the canopy
    !! \param units meters (m)
    real(rk)       ::    hgtref

    !> \brief Reference soil type
    !> \details Soil type
    integer        ::    sotypref

    !> \brief Reference surface pressure
    !> \details Surface pressure
    !! \param units hPa
    real(rk)       ::    pressfcref

    !> \brief Reference downward shortwave radiation
    !> \details Instantaneous downward shortwave radiation
    !! \param units W/m²
    real(rk)       ::    dswrfref

    !> \brief Reference surface sensible heat flux
    !> \details Instantaneous surface sensible heat net flux
    !! \param units W/m²
    real(rk)       ::    shtflref

    !> \brief Reference surface temperature
    !> \details Surface temperature
    !! \param units K
    real(rk)       ::    tmpsfcref

    !> \brief Reference 2-meter temperature
    !> \details 2-meter temperature
    !! \param units K
    real(rk)       ::    tmp2mref

    !> \brief Reference 2-meter specific humidity
    !> \details 2-meter specific humidity
    !! \param units kg/kg
    real(rk)       ::    spfh2mref

    !> \brief Reference planetary boundary layer height
    !> \details Height of planetary boundary layer
    !! \param units meters (m)
    real(rk)       ::    hpblref

    !> \brief Reference precipitation rate
    !> \details Mass precipitation rate
    !! \param units kg/m²/s
    real(rk)       ::    prate_averef

    !> \brief Reference soil moisture layer 1
    !> \details Volumetric soil moisture layer 1
    !! \param units m³/m³
    real(rk)       ::    soilw1ref

    !> \brief Reference soil moisture layer 2
    !> \details Volumetric soil moisture layer 2
    !! \param units m³/m³
    real(rk)       ::    soilw2ref

    !> \brief Reference soil moisture layer 3
    !> \details Volumetric soil moisture layer 3
    !! \param units m³/m³
    real(rk)       ::    soilw3ref

    !> \brief Reference soil moisture layer 4
    !> \details Volumetric soil moisture layer 4
    !! \param units m³/m³
    real(rk)       ::    soilw4ref

    !> \brief Reference wilting point
    !> \details Wilting point
    !! \param units proportion
    real(rk)       ::    wiltref

    !> \brief Reference ozone W126 values
    !> \details Ozone W126 values
    !! \param units ppm-hours
    real(rk)       ::    ozone_w126ref

    !> \brief Reference soil temperature level 1
    !> \details Soil temperature level 1
    !! \param units K
    real(rk)       ::    soilt1ref

    !> \brief Reference soil temperature level 2
    !> \details Soil temperature level 2
    !! \param units K
    real(rk)       ::    soilt2ref

    !> \brief Reference soil temperature level 3
    !> \details Soil temperature level 3
    !! \param units K
    real(rk)       ::    soilt3ref

    !> \brief Reference soil temperature level 4
    !> \details Soil temperature level 4
    !! \param units K
    real(rk)       ::    soilt4ref

    !> \brief Reference first model layer air temperature
    !> \details 1st model layer air temperature above ground
    !! \param units K
    real(rk)       ::    tmp_hyblev1ref

    !> \brief Reference average ground snow cover
    !> \details Average percent ground snow cover
    !! \param units percent (%)
    real(rk)       ::    snowc_averef

    !> \brief Reference average ground or water ice cover
    !> \details Average fraction ground or water ice cover
    !! \param units dimensionless fraction
    real(rk)       ::    icec_averef

!    real(rk)       ::    lev01ref, lev02ref, lev03ref, lev04ref, lev05ref, & !Input canopy profile levels
!                         lev06ref, lev07ref, lev08ref, lev09ref, lev10ref, &
!                         lev11ref, lev12ref, lev13ref, lev14ref
!    real(rk)       ::    pavd01ref, pavd02ref, pavd03ref, pavd04ref, pavd05ref, & !Input canopy PAVD profile
!                         pavd06ref, pavd07ref, pavd08ref, pavd09ref, pavd10ref, &
!                         pavd11ref, pavd12ref, pavd13ref, pavd14ref
    !> \brief Reference plant area volume density array
    !> \details Plant area volume density
    !! \param units m²/m³
    real(rk), allocatable   :: pavdref ( : ), pavd_arr ( : )

    !> \brief Reference vertical levels array
    !> \details Reference vertical levels with 3D input data
    !! \param units meters (m)
    real(rk), allocatable   :: levref ( : ), lev_arr  ( : )

    !> \}

END MODULE canopy_canmet_mod
