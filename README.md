<h1>
  <a href="https://github.com/noaa-oar-arl/canopy-app">
    <img src="docs/canopy-app-logo_no-bg.png" alt="canopy-app logo" height="125" valign="bottom">
  </a>
</h1>

[![License](https://img.shields.io/github/license/noaa-oar-arl/canopy-app.svg)](https://github.com/noaa-oar-arl/canopy-app/blob/main/LICENSE)
[![CI status](https://github.com/noaa-oar-arl/canopy-app/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/noaa-oar-arl/canopy-app/actions/workflows/ci.yml)

Repository for low-level, stand-alone/column canopy parameterizations for testing and application to gridded atmospheric composition/air quality models.

Authors: Patrick Campbell, Zachary Moon, and Wei-Ting Hung

## Getting Started

### Build

Canopy-App requires NetCDF-Fortran Libraries (i.e., `-lnetcdf -lnetcdff`) when using the 2D NetCDF I/O Option (i.e., `infmt_opt=0`).
See [the included Makefile](./src/Makefile), which detects NetCDF using `nf-config`, for an example (on GMU Hopper, you can use the `netcdf-c/4.7.4-vh` and `netcdf-fortran/4.5.3-ff` modules).

Compilation options can be controlled with environment variables:
 - `DEBUG=0` (off; default) or `DEBUG=1` (on)
 - `NETCDF=0` (off) or `NETCDF=1` (on; default)

Example:
```
DEBUG=1 NETCDF=1 make -C src
```

### Modify settings

If necessary, modify [the settings](#table-3-current-user-namelist-options) in the Fortran namelist file [`input/namelist.canopy`](./input/namelist.canopy),
which is read at runtime.

### Run

```
./canopy
```

## Components

Current Canopy-App components:

1.  In-Canopy Winds and Wind Adjustment Factor (WAF) for wildfire spread and air quality applications.  Based on Massman et al. (2017).

    Namelist Option : `ifcanwind` and/or `ifcanwaf` Output Variables: `canwind` (m s-1) `waf` (fraction)

    - `canopy_wind_mod.F90`
    - `canopy_waf_mod.F90`

2.  In-Canopy vertical diffusion (i.e., eddy diffusivities used to scale resolved model layer 1 diffusion).  Based on Massman et al. (2017) and Makar et al. (2017).

    Namelist Option : `ifcaneddy`  Output Variables:  `kz` (m2 s-1)

    - `canopy_eddyx_mod.F90`

3.  In-Canopy photolysis attenuation (i.e., used to scale resolved model layer 1 photolysis).  Based on Massman et al. (2017) and Markar et al. (2017).  

    Namelist Option : `ifcanphot`  Output Variables: `rjcf` (fraction)

    - `canopy_phot_mod.F90`

4.  In-Canopy leaf-level biogenic emissions (kg m-3 s-1). Based on MEGANv2 and v3 (Guenther et al., 2012), and using both Clifton et al. (2021) and Silva et al. (2020) parameterizations.

    - Note the emissions here are at leaf-level and the units are in per m3 (in each canopy layer volume using the LAD/biomass distribution) for the respective vegetation type in each grid cell/point. This is different then MEGANv2 or v3, as such models approximate combined activity factors per canopy level, sum them weighted to a given biomass distribution, and use total LAI to calculate the "big-leaf" 2D flux of biogenic emissions to the overlying atmosphere.  Thus, to get 2D total flux of biogenic emissions per m2 from Canopy-App, the explicit leaf-level emissions profile must be integrated across the individual canopy layer depths/resolutions (i.e., "modres", see [Table 3](#table-3-current-user-namelist-options)).  Since the modres is constant in Canopy-App, the layers can be directly summed and then multiplied by modres to get kg m-2 s-1. This sum neglects any effects of time integrated losses and/or chemistry that would reduce total biogenic emissions flux from the canopy.  When fractional vtypes (i.e., land use) are used, the summed layers can be multiplied by the fractional grid box areal coverage for each vegetation types in the grid cell.  However, for the current dominant vtype approach as input to Canopy-App and many UFS applications this multiplicative fraction = 1.

    Namelist Option : `ifcanbio`   Output Variables: see [Table 1](#table-1-canopy-app-biogenic-emissions-output-variables) below

    - `canopy_bioemi_mod.F90`

## Outputs

**Note for Biogenic emissions:** When `ifcanbio=.TRUE.`, output will include 3D canopy resolved biogenic emissions for the following species (based on Guenther et al., 2012), which have been mapped from Guenther et al. PFTs to input LU_OPT.

### Table 1. Canopy-App Biogenic Emissions Output Variables

| Variable Name    | Variable Description (Units: kg m-3 s-1)          |
| ---------------  | ------------------------------------------------- |
| emi_isop         | Isoprene                                          |
| emi_myrc         | Myrcene                                           |
| emi_sabi         | Sabinene                                          |
| emi_limo         | Limonene                                          |
| emi_care         | 3-Carene                                          |
| emi_ocim         | t-beta-Ocimene                                    |
| emi_bpin         | beta-Pinene                                       |
| emi_apin         | alpha-Pinene                                      |
| emi_mono         | Other Monoterpenes (34 compounds, Table 1 Guenther et al. (2012)               |
| emi_farn         | alpha-Farnesene                                   |
| emi_cary         | beta-Caryophyllene                                |
| emi_sesq         | Other Sesquiterpene (30 compounds, Table 1 Guenther et al. (2012)              |
| emi_mbol         | 232-MBO emissions                                 |
| emi_meth         | Methanol emissions                                |
| emi_acet         | Acetone emissions                                 |
| emi_co           | Carbon Monoxide emissions                         |
| emi_bvoc         | Bi-Directional VOC emissions (5 compounds, Table 1 Guenther et al. (2012)      |
| emi_svoc         | Stress VOC emissions (15 compounds, Table 1 Guenther et al. (2012)             |
| emi_ovoc         | Other VOC emissions (49 compounds, Table 1 Guenther et al. (2012)              |

**Current Canopy-App Output:** As discussed above, the current Canopy-App optional outputs includes 3D canopy winds (`canwind`), canopy vertical/eddy diffusivity values `kz`), biogenic emissions (see Table 1), and
canopy photolysis attenuation correction factors (`rjcf`).  Current 2D fields includes the Wind Adjustment Factor (`waf`).

Namelist Option : `file_out`  Prefix string (e.g., `'test'`) used to name output file (Output is 1D txt when using input 1D data (i.e., `infmt_opt=1`), or is 2D NetCDF output when 2D NetCDF input is used (i.e., `infmt_opt=0`)).

## Inputs and Settings

**Current Canopy-App Input:** Typical 1D or 2D (time=1,lat,lon) gridded atmospheric model input variables in 1st layer above canopy

Namelist Option : `file_vars`  Full name of input file (Supports either text or NetCDF format with following formats: `.txt`, `.nc`, `.ncf`, or `.nc4`)

- See example file inputs for variables and format (`gfs.t12z.20220701.sfcf000.canopy.txt` or `gfs.t12z.20220701.sfcf000.canopy.nc`).  Example surface met/land/soil inputs are based on NOAA's UFS-GFSv16 inputs initialized on July 01, 2022 @ 12 UTC (forecast at hour 000). Other external inputs for canopy related and other calculated variables are from numerous sources.  See [Table 2](#table-2-canopy-app-required-input-variables) below for more information.  **Note:** The example GFSv16 domain has been cut to the southeast U.S. region only in this example for size/time constraints here.
- Canopy-App assumes the NetCDF input files are in CF-Convention and test file is based on UFS-GFSv16; recommend using double or float for real variables.  Input data must be valid values.
- Canopy-App can also be run with a single point of 1D input data in a text file (e.g. `input_variables_point.txt`).

The Canopy-App input data in [Table 2](#table-2-canopy-app-required-input-variables) below is based around NOAA's UFS operational Global Forecast System Version 16 (GFSv16) gridded met data, and is supplemented with external canopy data (from numerous sources) and other external and calculated input variables.  

### Table 2. Canopy-App Required Input Variables

| **GFS /Met/Land/Soil Variables**    | **Variable Description and Units**                           |  **Example Data Sources/References (if necessary)**                                                |
| ----------------------------------  | ------------------------------------------------------------ |  ------------------------------------------------- 						      |
| lat                                 | Latitude  (degrees)                                          |  N/A                                               						      |
| lon                                 | Longitude (degrees; from 0-360)                              |  N/A                                               					              |
| time                                | Timestamp (days since YYYY-N-D 0:0:0) (NetCDF Only)          |  N/A                                               					              |  
| ugrd10m                             | U wind at reference height above canopy (m/s), e.g., 10 m    |  UFS NOAA/GFSv16 *(see below for downloading using AWS)                                 						      |
| vgrd10m                             | V wind at reference height above canopy  (m/s), e.g., 10 m   |  UFS NOAA/GFSv16                                   						      |
| vtype                               | Vegetation type (dimensionless), VIIRS or MODIS              |  UFS NOAA/GFSv16                                   						      |
| fricv                               | Friction velocity (m/s)                                      |  UFS NOAA/GFSv16                                   						      |
| sfcr                                | Total surface roughness length (m)                           |  UFS NOAA/GFSv16                                   						      |
| sotyp            			  | Soil type (dimensionless), STATSGO                           |  UFS NOAA/GFSv16                                   						      |
| pressfc          			  | Surface pressure (Pa)                                        |  UFS NOAA/GFSv16                                   						      |
| dswrf                               | Instantaneous downward shortwave radiation at surface (W/m2) |  UFS NOAA/GFSv16                                   						      |
| shtfl                               | Instantaneous sensible heat flux at surface (W/m2)           |  UFS NOAA/GFSv16                                   						      |
| tmpsfc                              | Surface temperature (K)                                      |  UFS NOAA/GFSv16                                  						      |
| tmp2m                               | 2-meter temperature (K)                                      |  UFS NOAA/GFSv16                                   						      |
| spfh2m                              | 2-meter specific humidity (kg/kg)                            |  UFS NOAA/GFSv16                                   					              |
| hpbl                                | Height of the planetary boundary layer (m)                   |  UFS NOAA/GFSv16                                   						      |
| prate_ave                           | Average mass precipitation rate (kg m-2 s-1)                 |  UFS NOAA/GFSv16                                   						      |
| **External Canopy Variables**       | **Variable Description and Units**                           |  **Data Source/Reference (if necessary)**         						      |
| fh                                  | Forest canopy height (m)                                     |  Fused GEDI/Landsat data. Data Period=2020. ([Potapov et al., 2020](https://doi.org/10.1016/j.rse.2020.112165))           |  
| clu                                 | Forest clumping index (dimensionless)                        |  GridingMachine/MODIS. Data Period=2001-2017 Climatology. ([Wei et al., 2019](https://doi.org/10.1016/j.rse.2019.111296)). Extended globally for high latitudes using methods described [here](https://gmuedu-my.sharepoint.com/:w:/g/personal/whung_gmu_edu/EdglXmW2kzBDtDj1xV0alGcB1Yo2I8hzdyWGVGB2YOTfgw). |
| lai                                 | Leaf area index (m2/m2)                                      |  VIIRS-NPP. Data Period=2018-2020 Climatology. ([Myneni 2018](https://doi.org/10.5067/VIIRS/VNP15A2H.001)). Extended globally for high latitudes using methods described [here](https://gmuedu-my.sharepoint.com/:w:/g/personal/whung_gmu_edu/EdglXmW2kzBDtDj1xV0alGcB1Yo2I8hzdyWGVGB2YOTfgw). |                |
| ffrac                               | Forest fraction (dimensionless)                              |  Based on [VIIRS GVF-NPP](https://www.star.nesdis.noaa.gov/jpss/gvf.php) and GriddingMachine/MODIS FFRAC/FVC from Terra.  Data Period=2020. ([DiMiceli et al., 2022](https://doi.org/10.5067/MODIS/MOD44B.061)). Extended globally for high latitudes using methods described [here](https://gmuedu-my.sharepoint.com/:w:/g/personal/whung_gmu_edu/EdglXmW2kzBDtDj1xV0alGcB1Yo2I8hzdyWGVGB2YOTfgw). |      |  
| **Other External Variables**        | **Variable Description and Units**                           |  **Data Source/Reference (if necessary)**                                                          |
| frp                                 | Total Fire Radiative Power (MW/grid cell area)               |  [NOAA/NESDIS GBBEPx](https://www.ospo.noaa.gov/Products/land/gbbepx/)                              |
| csz                                 | Cosine of the solar zenith angle (dimensionless)             |  [Based on Python Pysolar](https://pysolar.readthedocs.io/en/latest/)                               |
| mol                                 | Monin-Obukhov Length (m)                                     |  Externally calculated using GFS tmp2m, fricv, and shtfl.  ([Essa, 1999](https://inis.iaea.org/collection/NCLCollectionStore/_Public/37/118/37118528.pdf))   |
| href                                | Reference height above canopy (m) - 10 m                     |  Added as a constant reference height above surface (i.e., 10 m).  Can be taken from NL.           |

**More Information on Data Sources from [Table 2](#table-2-canopy-app-required-input-variables):**

**Downloading GFS Files from AWS:** NOAA's hourly global GFS, gridded (at ~13x13 km resolution) data may be downloaded publicly from the following Amazon Web Service (AWS) S3 location:
```
https://nacc-in-the-cloud.s3.amazonaws.com/inputs/YYYYMMDD/gfs.t12z.sfcfHHH.nc
```
Where HHH pertains to the hour of the 24-hr forecast (e.g., f000 is initialization).  
Example download command using wget:
```
wget --no-check-certificate --no-proxy https://nacc-in-the-cloud.s3.amazonaws.com/inputs/20230215/gfs.t12z.sfcf000.nc
```
Hourly gridded GFSv16 data is available on AWS from March 23, 2021 - Current Day.

**GriddingMachine:** GriddingMachine is open source database and software for Earth system modeling at global and regional scales.  Data is easily accessible in consistent formats for ease of downloading/processing.  All available datasets may be found at:  https://github.com/CliMA/GriddingMachine.jl. ([Wang et al., 2022](https://doi.org/10.1038/s41597-022-01346-x)).

**Downloading Example Canopy Files from AWS:** Example monthly, global gridded files containing all GFSv16 met/land/soil data from 2022 combined with external canopy and other external variables (regridded to GFSv16 13 km resolution) described above may also be downloaded via AWS S3 location:
```
https://nacc-in-the-cloud.s3.amazonaws.com/inputs/geo-files/canopy.2022MM01.testf000.global.nc
```

### Table 3. Current User Namelist Options

| Namelist Option  | Namelist Description and Units                                                     |
| ---------------  | ---------------------------------------------------------------------------------- |
| infmt_opt        | integer for choosing 1D text ( = 1)  or 2D NetCDF input file format (= 0, default) |
| nlat             | number of latitude cells (must match # of LAT in `file_vars` above)                |
| nlon             | number of longitude cells (must match # of LON in `file_vars` above)               |
| modlays          | number of model (below and above canopy) layers                                    |
| modres           | above and below canopy model vertical resolution (m)                               |
| ifcanwind        | logical canopy wind option (default = .FALSE.)                                     |
| ifcanwaf         | logical canopy WAF option (default = .FALSE.) **\*\***                             |
| ifcaneddy        | logical canopy eddy Kz option (default = .FALSE.)                                  |
| ifcanphot        | logical canopy photolysis option (default = .FALSE.)                               |
| ifcanbio         | logical canopy biogenic emissions option (default = .FALSE.)                       |
| href_opt         | integer for using href_set in namelist (= 0) or array from file (= 1); default = 0 |
| href_set         | user set real value of reference height above canopy associated with input wind speed (m) (only used if href_opt = 0)  |
| z0ghc            | ratio of ground roughness length to canopy top height (Massman et al., 2017)       |
| lambdars         | Value representing influence of roughness sublayer (Massman et al., 2017)          |
| dx_opt           | 0=Calculation of dx resolution/distance from lon; 1=user set dx grid resolution    |
| dx_set           | user set real value of grid resolution (m) only if dx_opt=1                        |
| flameh_opt       | 0=Calculation of flame height from FRP (Byram, 1959); 1=user set flameh; 2=FRP calculation where available (active fires), elsewhere user set flameh; 3=FlameH override, i.e., only uses fraction of canopy height (flameh_set must be <=1.0) as a surrogate for flameh; 4=FRP calculation where available (active fires) and FlameH override elsewhere (same as option 3); 5=FRP/intensity dependent (i.e., sub-canopy vs. crown fires) calculation where available (active fires) and FlameH override elsewhere (same as option 3). If option 5 is used and crowning is calculated, then the total flame height (i.e., top of canopy=FCH) is used instead of 1/2 flame height.    |
| flameh_set       | user set real value of flame height (m) if flameh_opt = 1 or 2, or flameh = fraction of canopy height (<=1.0), i.e., flameh override, if flameh_opt = 3, 4, or 5                       |
| frp_fac          | user set real value of tuning factor applied to FRP in calculation of flame height (Default = 1.0). Used only if flameh_opt = 0, 2, 4, or 5.                      |
| pai_opt          | integer (0=PAI fixed from Katul et al. 2004 veg types-->default; 1=PAI Massman et al. 2017 Eq. 19 calc; 2=PAI from model LAI+WAI; 3=user set PAI value) |
| pai_set          | user set real value of PAI (default=4.0; only used if pai_opt=3)                   |
| lu_opt           | integer for input model land use type (0=VIIRS 17 Cat (default) or 1=MODIS-IGBP 20 Cat (valid LU types 1-10 and 12); input mapped to Massman et al.)   |
| z0_opt           | integer (0=use model input or 1 = vegtype dependent z0 for first estimate)         |
| rsl_opt          | user set option to include more explicit stability and Roughness SubLayer (RSL) effects in calculation of wind speed at canopy top (Uc) from reference height  0 = off; 1 = on |
| bio_cce          | user set real value of MEGAN biogenic emissions "canopy environment coefficient" used to tune emissions to model inputs/calculations (default = 0.21 based on Silva et al. 2020)   |
| lai_thresh       | user set real value of LAI threshold for contiguous canopy (m2/m2)                 |
| frt_thresh       | user set real value of forest fraction threshold for contiguous canopy             |
| fch_thresh       | user set real value of canopy height  threshold for contiguous canopy (m)          |

**\*\*** If modres >> flameh then some error in WAF calculation will be incurred.  Suggestion is to use relative fine modres (at least <= 0.5 m) compared to average flame heights (e.g., ~ 1.0 m) if WAF is required.

**Note:** Canopy is parameterized by foliage distribution shape functions and parameters for different vegetation types.

- `canopy_profile_mod.F90`

## References

*Further references contained within the source code.*

- Clifton, O. E., Patton, E. G., Wang, S., Barth, M., Orlando, J., & Schwantes, R. H. (2022). Large eddy simulation for investigating coupled forest canopy and turbulence influences on atmospheric chemistry. Journal of Advances in Modeling Earth Systems, 14, e2022MS003078. https://doi. org/10.1029/2022MS003078

- Guenther, A. B., Jiang, X., Heald, C. L., Sakulyanontvittaya, T., Duhl, T., Emmons, L. K., and Wang, X.: The Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling biogenic emissions, Geosci. Model Dev., 5, 1471–1492, https://doi.org/10.5194/gmd-5-1471-2012, 2012.

- Katul, G.G., Mahrt, L., Poggi, D., and Sanz, C. (2004). One- and two-equation models for canopy turbulence. Boundary-Layer Meteorol. 113: 81–109. https://doi.org/10.1023/B:BOUN.0000037333.48760.e5

- Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone. Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms1524

- Massman, W. J., J.M. Forthofer, and M.A. Finney. (2017). An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior. Canadian Journal of Forest Research. 47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354

- Silva, S. J., Heald, C. L., and Guenther, A. B.: Development of a reduced-complexity plant canopy physics surrogate model for use in chemical transport models: a case study with GEOS-Chem v12.3.0, Geosci. Model Dev., 13, 2569–2585, https://doi.org/10.5194/gmd-13-2569-2020, 2020.

## Development

After cloning the repository,
set up the pre-commit hooks by invoking
```
pre-commit install --install-hooks
```
within the repository (after [installing `pre-commit`](https://pre-commit.com/#installation)).
This only has to be done once.
The [current configuration](./.pre-commit-config.yaml) applies
[`findent -i4`](https://www.ratrabbit.nl/ratrabbit/findent/) to fix indentation
and strips trailing whitespace.
Using pre-commit saves you time/energy and reduces diffs.
