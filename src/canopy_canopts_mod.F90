!> \file canopy_canopts_mod.F90
!> \brief Canopy model configuration options and user-settable parameters
!> \details This module contains all user-configurable options and parameters
!>          for the canopy model, including vegetation types, emission options,
!>          dry deposition settings, and model physics options.
!> \author P. C. Campbell
!> \date 03 Oct 2022

!> \defgroup canopy_options Canopy Model Options
!> \brief Configuration options and user-settable parameters for canopy model

MODULE canopy_canopts_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Option Variable Descriptions
! Purpose:  Contains canopy option variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk
    IMPLICIT NONE

!! .... defines canopy options (read from user namelist)
    !> \ingroup canopy_options
    !> \{

    !> \brief Input file format option
    !> \details Integer for choosing 1D or 2D input file format
    !! - 0 = 2D format (default)
    !! - 1 = 1D format
    integer             ::    infmt_opt

    !> \brief 3D variable input option
    !> \details Integer for choosing if 3D variables will be read from file
    !! - 0 = off (default)
    !! - 1 = on
    integer             ::    var3d_opt

    !> \brief Number of 3D levels in input file
    !> \details Integer for number of 3D levels in input file, if var3d_opt=1
    !! \param default 14 levels
    integer             ::    var3d_set

    !> \brief 3D GEDI PAVD profiles option
    !> \details Integer for turning on 3D GEDI PAVD profiles read from input file
    !! - 0 = off (default)
    !! - 1 = on
    integer             ::    pavd_opt

    !> \brief GEDI PAVD latitude threshold
    !> \details Real value for latitude +/- threshold when 3D GEDI PAVD profiles read from input file
    !! \param default 52.0 degrees
    real(rk)            ::    pavd_set

    !> \brief Reference height option
    !> \details Integer for using set href in namelist (=0) or array from file(=1)
    !! - 0 = use namelist value (default)
    !! - 1 = read from file
    integer             ::    href_opt

    !> \brief Set reference height above canopy
    !> \details Reference height above canopy at 10 m
    !! \param units meters (m)
    !! \param default 10.0 m
    real(rk)            ::    href_set

    !> \brief Canopy wind option
    !> \details Logical canopy wind option
    !! \param default .FALSE.
    logical             ::    ifcanwind

    !> \brief Canopy WAF option
    !> \details Logical canopy WAF (Wind Attenuation Factor) option
    !! \param default .FALSE.
    logical             ::    ifcanwaf

    !> \brief Canopy eddy Kz option
    !> \details Logical canopy eddy diffusivity Kz option
    !! \param default .FALSE.
    logical             ::    ifcaneddy

    !> \brief Canopy photolysis attenuation option
    !> \details Logical canopy photolysis attenuation option
    !! \param default .FALSE.
    logical             ::    ifcanphot

    !> \brief Canopy biogenic emissions option
    !> \details Logical canopy biogenic emissions option
    !! \param default .FALSE.
    logical             ::    ifcanbio

    !> \brief Canopy gas dry deposition option
    !> \details Logical canopy gas dry deposition option
    !! \param default .FALSE.
    logical             ::    ifcanddepgas
    !> \brief PAI calculation option
    !> \details Integer for PAI (Plant Area Index) values used or calculated
    !! - 0 = calculated from model (default)
    !! - 1 = use set value
    integer             ::    pai_opt

    !> \brief Set PAI value
    !> \details Real value for PAI set values used when pai_opt=1
    !! \param units m²/m²
    !! \param default 4.0
    real(rk)            ::    pai_set

    !> \brief Land use mapping option
    !> \details Integer for LU type from model mapped to Massman et al.
    !! - 0 = VIIRS (default)
    !! - 1 = other mapping
    integer             ::    lu_opt

    !> \brief Surface roughness option
    !> \details Integer for setting first estimate of z0
    !! - 0 = Z0_MOD (default)
    !! - 1 = other method
    integer             ::    z0_opt

    !> \brief Flame height calculation option
    !> \details Integer for flame height values used or calculated
    !! - 0 = calculated (default)
    !! - 1 = use set value
    integer             ::    flameh_opt

    !> \brief FRP to flame height relationship
    !> \details Integer for FRP (Fire Radiative Power) to flame height relationships used
    !! - 0 = default relationship
    !! - 1 = alternative relationship
    integer             ::    flameh_cal

    !> \brief User set flame height
    !> \details User Set Flame Height value
    !! \param units meters (m)
    real(rk)            ::    flameh_set

    !> \brief FRP tuning factor
    !> \details FRP tuning factor for flame height calculation
    !! \param default 1.0
    real(rk)            ::    frp_fac

    !> \brief Grid resolution option
    !> \details Integer for dx resolution values used or calculated
    !! - 0 = calculated (default)
    !! - 1 = use set value
    integer             ::    dx_opt

    !> \brief User set grid cell resolution
    !> \details User Set Grid Cell Resolution
    !! \param units meters (m)
    real(rk)            ::    dx_set

    !> \brief LAI threshold for canopy conditions
    !> \details User set grid cell LAI threshold to apply canopy conditions
    !! \param units m²/m²
    real(rk)            ::    lai_thresh

    !> \brief Canopy fraction threshold
    !> \details User set grid cell canopy fraction threshold to apply canopy conditions
    !! \param units dimensionless fraction
    real(rk)            ::    cf_thresh

    !> \brief Canopy height threshold
    !> \details User set grid cell canopy height threshold to apply canopy conditions
    !! \param units meters (m)
    real(rk)            ::    ch_thresh

    !> \brief Roughness sublayer option
    !> \details RSL option used in model from Rosenzweig et al. 2021
    !! - 0 = off (default)
    !! - 1 = on
    integer             ::    rsl_opt

    !> \brief Ground roughness to canopy height ratio
    !> \details Ratio of ground roughness length to canopy top height
    !! \param units dimensionless
    real(rk)            ::    z0ghc

    !> \brief Roughness sublayer influence parameter
    !> \details Value representing influence of roughness sublayer
    !! \param units dimensionless
    real(rk)            ::    lambdars
    !> \brief MEGAN biogenic emission canopy environment coefficient
    !> \details MEGAN biogenic emission canopy environment coefficient
    !! \param units dimensionless
    real(rk)            ::    bio_cce

    !> \brief Biogenic species output option
    !> \details Set default integer for species output option
    !! - 0 = all species (default)
    !! - 1 = selected species
    integer             ::    biospec_opt

    !> \brief MEGAN vertical integration option
    !> \details MEGAN vertical integration of emissions option
    !! - 0 = off (default)
    !! - 1 = on
    integer             ::    biovert_opt

    !> \brief Canopy data source option
    !> \details Set default integer for canopy option from GEDI or user
    !! - 0 = GEDI data (default)
    !! - 1 = user-defined values
    integer             ::    can_opt

    !> \brief Default canopy vegetation height
    !> \details Set default value for canopy vegtype heights used in model
    !! \param units meters (m)
    !! \param default 1.0 m
    real(rk)            ::    can_chset

    !> \brief Default canopy vegetation fraction
    !> \details Set default value for canopy vegfrac used in model
    !! \param units dimensionless fraction
    !! \param default 0.5
    real(rk)            ::    can_cfset

    !> \brief Default canopy LAI
    !> \details Set default value for canopy LAI used in model
    !! \param units m²/m²
    !! \param default 0.1
    real(rk)            ::    can_laiset

    !> \brief Shrubs/savanna/grassland vegetation option
    !> \details Set default integer for shrubs/savanna/grassland vegtype option from GEDI or user
    !! - 0 = GEDI data (default)
    !! - 1 = user-defined values
    integer             ::    ssg_opt

    !> \brief Default shrub/savanna/grassland height
    !> \details Set default value for shrubs/savanna/grassland vegtype heights used in model
    !! \param units meters (m)
    !! \param default 1.0 m
    real(rk)            ::    ssg_chset

    !> \brief Default shrub/savanna/grassland vegetation fraction
    !> \details Set default value for shrubs/savanna/grassland vegfrac used in model
    !! \param units dimensionless fraction
    !! \param default 0.5
    real(rk)            ::    ssg_cfset
    !> \brief Default shrub/savanna/grassland LAI
    !> \details Set default value for shrubs/savanna/grassland LAI used in model
    !! \param units m²/m²
    !! \param default 0.1
    real(rk)            ::    ssg_laiset

    !> \brief Crop vegetation option
    !> \details Set default integer for crop vegtype option from GEDI or user
    !! - 0 = GEDI data (default)
    !! - 1 = user-defined values
    integer             ::    crop_opt

    !> \brief Default crop height
    !> \details Set default value for crop vegtype heights used in model
    !! \param units meters (m)
    !! \param default 3.0 m
    real(rk)            ::    crop_chset

    !> \brief Default crop vegetation fraction
    !> \details Set default value for crop vegfrac used in model
    !! \param units dimensionless fraction
    !! \param default 0.5
    real(rk)            ::    crop_cfset

    !> \brief Default crop LAI
    !> \details Set default value for crop LAI used in model
    !! \param units m²/m²
    !! \param default 0.1
    real(rk)            ::    crop_laiset

    !> \brief CO₂ inhibition option for isoprene
    !> \details Set default integer for CO₂ inhibition option for biogenic isoprene emissions
    !! - 0 = Possell & Hewitt (2011) (default)
    !! - 1 = alternative method
    integer             ::    co2_opt

    !> \brief Atmospheric CO₂ concentration
    !> \details Set default value for atmospheric CO₂ concentration for co2_opt
    !! \param units ppmv
    !! \param default 400.0 ppmv
    real(rk)            ::    co2_set

    !> \brief Leaf age factor option
    !> \details Set default for Leaf Age factor option for BVOCs
    !! - 1 = GAMMA_LEAFAGE = 1 (default)
    integer             ::    leafage_opt

    !> \brief LAI timestep
    !> \details Set default timestep for LAI input
    !! \param units seconds
    !! \param default 86400 (daily), 2592000 for monthly
    integer             ::    lai_tstep

    !> \brief Canopy loss ratios option
    !> \details Set default integer for turning on canopy loss ratios for adjusting top of canopy net emissions
    !! - 0 = off (default)
    !! - 1 = use lifetime-based loss
    !! - 2 = use constant loss factor
    integer             ::    loss_opt

    !> \brief Above-canopy BVOC lifetime
    !> \details Set default value for above-canopy BVOC lifetime used with loss_opt=1
    !! \param units seconds (s)
    !! \param default 3600 s
    real(rk)            ::    lifetime

    !> \brief Constant canopy loss factor
    !> \details Set default value for constant canopy loss factor applied used with loss_opt=2
    !! \param units dimensionless
    !! \param default 0.96
    real(rk)            ::    loss_set

    !> \brief Loss factor application index
    !> \details Set default integer for applying canopy loss factor to all species (=0) or only specific biogenics specie indices (> 0)
    integer             ::    loss_ind

    !> \brief Historical averaging option
    !> \details Set default integer for using historically averaged leaf temp and PAR for biogenic emissions
    !! - 0 = off (default)
    !! - 1 = on
    integer             ::    hist_opt

    !> \brief Soil moisture factor option
    !> \details Set default for Soil Moisture factor option for BVOCs
    !! - 1 = GAMMA_SOIM = 1 (default)
    integer             ::    soim_opt
    !> \brief Soil layer 1 depth
    !> \details User set real value of depth of soil layer 1 centerpoint
    !! \param units centimeters (cm)
    !! \param default 5 cm (based on Noah/Noah-MP)
    real(rk)            ::    soild1

    !> \brief Soil layer 2 depth
    !> \details User set real value of depth of soil layer 2 centerpoint
    !! \param units centimeters (cm)
    !! \param default 25 cm (based on Noah/Noah-MP)
    real(rk)            ::    soild2

    !> \brief Soil layer 3 depth
    !> \details User set real value of depth of soil layer 3 centerpoint
    !! \param units centimeters (cm)
    !! \param default 70 cm (based on Noah/Noah-MP)
    real(rk)            ::    soild3

    !> \brief Soil layer 4 depth
    !> \details User set real value of depth of soil layer 4 centerpoint
    !! \param units centimeters (cm)
    !! \param default 150 cm (based on Noah/Noah-MP)
    real(rk)            ::    soild4

    !> \brief Air quality stress index option
    !> \details Set default integer for air quality stress index for gamma_aq in biogenic emissions
    integer             ::    aq_opt

    !> \brief Ozone W126 value
    !> \details Set default value for constant ozone W126 value
    !! \param units ppm-hours
    real(rk)            ::    w126_set

    !> \brief High temperature stress index option
    !> \details Set default integer for high temperature stress index for gamma_ht in biogenic emissions
    integer             ::    ht_opt

    !> \brief Low temperature stress index option
    !> \details Set default integer for low temperature stress index for gamma_lt in biogenic emissions
    integer             ::    lt_opt

    !> \brief High wind stress index option
    !> \details Set default integer for high wind stress index for gamma_hw in biogenic emissions
    integer             ::    hw_opt

    !> \brief Dry deposition species output option
    !> \details Set default integer for species output option
    !! - 0 = all species (default)
    !! - 1 = selected species
    integer             ::    ddepspecgas_opt

    !> \brief Chemical mechanism option
    !> \details Set default integer value to select chemical mechanism
    !! - 0 = RACM2 (default)
    integer             ::    chemmechgas_opt

    !> \brief Total chemical mechanism species
    !> \details Set default integer value to select chemical mechanism gas species list including transported
    !! \param default 31 for RACM2
    integer             ::    chemmechgas_tot

    !> \brief Soil category option
    !> \details Set default integer value to select soil category option
    !! - 0 = STATSGO/FAO (default)
    integer             ::    soilcat_opt

    !> \brief First hybrid model layer height
    !> \details Set default approximate input height of 1st hybrid model layer above ground (used for temp lapse rate approximation)
    !! \param units meters (m)
    !! \param default 20.0 meters
    real(rk)            ::    hyblev1

    !> \brief Snow cover threshold
    !> \details Set default snow cover percent at grid/point, above which ground surface is treated as dominant snow
    !! \param units percent (%)
    !! \param default 50%
    real(rk)            ::    snowc_set

    !> \brief Ice cover threshold
    !> \details Set default ice cover percent at grid/point, above which ground or water surface is treated as dominant ice
    !! \param units percent (%)
    !! \param default 50%
    real(rk)            ::    icec_set

    !> \brief Reaction probability for building surfaces
    !> \details Set default reaction probability for gas dry deposition to different building surfaces
    !! \param units dimensionless
    !! \param default 5.0D-5
    real(rk)            ::    gamma_set

    !> \brief Minimum aerodynamic resistance
    !> \details Set default minimum aerodynamic resistance
    !! \param units s/m
    !! \param default 10 s/m
    real(rk)            ::    Ramin_set

    !> \}
END MODULE canopy_canopts_mod
