# canopy-app
Repository for low-level, stand-alone/column canopy parameterizations for testing and application to gridded atmospheric composition/air quality models.

Author(s):

Patrick Campbell, Zachary Moon, and Wei-Ting Hung

Build canopy model:

Canopy-App requires NetCDF-Fortran Libraries (i.e., -lnetcdf -lnetcdff) when using the 2D NetCDF I/O Option (i.e., infmt_opt=0).
See included Makefile for example (currently developed using netcdf-c/4.7.4-vh and netcdf-fortran/4.5.3-ff modules).

Compilation Options with or without Debug or NetCDF options:
 - `DEBUG  = 0(off; default) or DEBUG  =1(on)`
 - `NETCDF = 0(off) or          NETCDF =1(on; default)`

Example: `DEBUG=1 NETCDF=1 make`

Compile, edit namelist, and run canopy model:
- `make`
- `namelist.canopy`
- `./canopy`

Canopy is parameterized by foliage distribution shape functions and parameters for different vegetation types.

- `canopy_parm_mod.F90`

Current Canopy-App components:

1.  In-Canopy Winds and Wind Adjustment Factor (WAF) for wildfire spread and air quality applications.

    Namelist Option : `ifcanwind`

    - `canopy_wind_mod.F90`
    - `canopy_waf_mod.F90`

2.  In-Canopy vertical diffusion (i.e., eddy diffusivity adjustment) for air quality applications.

    Namelist Option : `ifcaneddy`

    - `canopy_eddyx_mod.F90`

3.  In-Canopy vertical photolysis attenuation for air quality applications.

    Namelist Option : `ifcanphot`

    - `canopy_phot_mod.F90`

4.  In-Canopy vertical biogenic emissions (kg m-2 s-1) based on MEGANv2 and v3 and using Clifton et al. (2021) and Silva et al. (2020) parameterizations for air quality applications.

    Namelist Option : `ifcanbio`

    - `canopy_bioemi_mod.F90`

    **Current Canopy-App Input:** Typical 1D or 2D (time=1,lat,lon) gridded atmospheric model input variables in 1st layer above canopy

    Namelist Option : `file_vars`  Full name of input file (Supports either text or NetCDF format with following formats:
                                                            `.txt`, `.nc`, `.ncf`, or `.nc4`)

- See example file inputs for variables and format (`gfs.t12z.20220701.sfcf000.txt` or `gfs.t12z.20220701.sfcf000.nc`).  Example input is based on NOAA's UFS-GFSv16 inputs initialized on July 01, 2022 @ 12 UTC (forecast at hour 000)
- Canopy-App assumes the NetCDF input files are in CF-Convention and test file is based on UFS-GFSv16; recommend using double or float for real variables.
- Canopy-App can also be run with a single point of 1D input data in a text file (e.g. `input_variables_point.txt`).

    **Current Canopy-App Output:** Outputs 3D canopy winds, canopy vertical/eddy diffusivity values, and
    canopy photolysis attenuation correction factors, and 2D Wind Adjustment Factor (WAF).

    Namelist Option : `file_out`  Prefix string (e.g., 'test') used to name output file (Output is 1D txt when using input 1D data (i.e., infmt_opt=1), or is 2D NetCDF output when 2D NetCDF input is used (i.e., infmt_opt=0)).


    **Table 1. Canopy-App Required Input Variables**

    | Variable Name    | Variable Description and Units                    |
    | ---------------  | ------------------------------------------------- |  
    | lat              | Latitude  (degrees)                               |
    | lon              | Longitude (degrees; from 0-360)                   |
    | time             | Timestamp (days since YYYY-N-D 0:0:0) (NetCDF Only) |
    | fh               | Forest canopy height (m)                          |
    | href             | Reference height above canopy (m) - 10 m          |
    | ugrd10m          | U wind at HREF (m/s), e.g., 10 m                  |
    | vgrd10m          | V wind at HREF (m/s), e.g., 10 m                  |
    | clu              | Forest clumping index (dimensionless)             |
    | lai              | Leaf area index (m2/m2)                           |
    | vtype            | Vegetation type (dimensionless), VIIRS or MODIS   |
    | ffrac            | Forest fraction (dimensionless)                   |
    | fricv            | Friction velocity (m/s)                           |
    | csz              | Cosine of the solar zenith angle (dimensionless)  |
    | sfcr             | Total surface roughness length (m)                |
    | mol              | Monin-Obukhov Length (m)                          |
    | frp              | Total Fire Radiative Power (MW/grid cell area)    |
    | sotyp            | Soil type (dimensionless), STATSGO                |
    | pressfc          | Surface pressure (Pa)                             |
    | dswrf            | Instantaneous downward shortwave radiation at surface (W/m2) |
    | shtfl            | Instantaneous sensible heat flux at surface (W/m2)           |
    | tmpsfc           | Surface temperature (K)                                      |
    | tmp2m            | 2-meter temperature (K)                                      |
    | spfh2m           | 2-meter specific humidity (kg/kg)                            |
    | hpbl             | Height of the planetary boundary layer (m)                   |
    | prate_ave        | Average mass precipitation rate (kg m-2 s-1)                 |


    **Table 2. Current User Namelist Options**

    | Namelist Option  | Namelist Description and Units                                                     |
    | ---------------  | ---------------------------------------------------------------------------------- |
    | infmt_opt        | integer for choosing 1D text ( = 1)  or 2D NetCDF input file format (= 0, default) |
    | nlat             | number of latitude cells (must match # of LAT in `file_vars` above)                |
    | nlon             | number of longitude cells (must match # of LON in `file_vars`above)                |
    | modlays          | number of model (below and above canopy) layers                                    |
    | modres           | above and below canopy model vertical resolution (m)                               |
    | ifcanwind        | logical canopy wind option (default = .FALSE.)                                     |
    | ifcanwaf         | logical canopy WAF option (default = .FALSE.)**                                    |
    | ifcaneddy        | logical canopy eddy Kz option (default = .FALSE.)                                  |
    | ifcanphot        | logical canopy photolysis option (default = .FALSE.)                               |
    | ifcanbio         | logical canopy biogenic emissions option (default = .FALSE.)                       |
    | href_opt         | integer for using href_set in namelist (= 0) or array from file (= 1); default = 0 |
    | href_set         | user set real value of reference height above canopy associated with input wind speed (m) (only used if href_opt = 0)  |
    | z0ghc            | ratio of ground roughness length to canopy top height (Massman et al., 2017)       |
    | lambdars         | Value representing influence of roughness sublayer (Massman et al., 2017)          |
    | dx_opt           | 0=Calculation of dx resolution/distance from lon; 1=user set dx grid resolution    |
    | dx_set           | user set real value of grid resolution (m) only if dx_opt=1                        |
    | flameh_opt       | 0=Calculation of flame height from FRP (Byram, 1959); 1=user set flameh; 2=FRP calculation where available (active fires), elsewhere user set flameh    |
    | flameh_set       | user set real value of flame height (m) only if flame_opt=1                        |
    | pai_opt          | integer (0=PAI fixed from Katul et al. 2004 veg types-->default; 1=PAI Massman et al. 2017 Eq. 19 calc; 2=PAI from model LAI+WAI; 3=user set PAI value) |
    | pai_set          | user set real value of PAI (default=4.0; only used if pai_opt=3)                   |
    | lu_opt           | integer for input model land use type (0=VIIRS 17 Cat (default) or 1=MODIS-IGBP 20 Cat (valid LU types 1-10 and 12); input mapped to Massman et al.)   |
    | z0_opt           | integer (0=use model input or 1 = vegtype dependent z0 for first estimate)         |
    | lai_thresh       | user set real value of LAI threshold for contiguous canopy (m2/m2)                 |
    | frt_thresh       | user set real value of forest fraction threshold for contiguous canopy             |
    | fch_thresh       | user set real value of canopy height  threshold for contiguous canopy (m)          |
    | rsl_opt          | user set option to include more explicit stability and Roughness SubLayer (RSL) effects in calculation of wind speed at canopy top (Uc) from reference height  0 = off; 1 = on |

**Note:  If modres >> flameh then some error in WAF calculation will be incurred.  Suggestion is to use relative fine modres (at least <= 0.5 m) compared to average flame heights (e.g., ~ 1.0 m) if WAF is required.


Main Citations (further references contained within):

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
