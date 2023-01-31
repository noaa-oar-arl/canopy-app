# canopy-app
Repository for low-level, stand-alone/column canopy parameterizations for testing and application to gridded atmospheric composition/air quality models.

Author(s):

Patrick Campbell, Zachary Moon, and Wei-Ting Hung

Build canopy model:

Canopy-App requires NetCDF-Fortran Libraries (i.e., -lnetcdf -lnetcdff) for NetCDF I/O Option. See included Makefile for example.

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

    **Current Canopy-App Input:** Typical 1D or 2D (time=1,lat,lon) gridded atmospheric model input variables in 1st layer above canopy

    Namelist Option : `file_vars`  Full name of input file (Supports either text or NetCDF format with following formats:
                                                            `.txt`, `.nc`, `.ncf`, or `.nc4`)

- See example file inputs for variables and format (`input_variables.txt`, `input_variables_1D.nc`, or `input_variables_2D.nc`).
- Canopy-App assumes the NetCDF input files are in CF-Convention; recommend using double or float for real variables.
- Canopy-App can also be run with a single point of 1D input data in a text file (e.g. `input_variables_point.txt`).

    **Current Canopy-App Output:** Outputs 3D canopy winds, canopy vertical/eddy diffusivity values, and
    canopy photolysis attenuation correction factors, and 2D Wind Adjustment Factor (WAF).

    Namelist Option : `file_out`  Prefix string (e.g., 'test') used for output file name (Currently 1D txt output with input 1D data, or 2D NetCDF output when 2D NCF input is used, i.e., infmt_opt=0).


    **Table 1. Canopy-App Required Input Variables**

    | Variable Name    | Variable Description and Units                    |
    | ---------------  | ------------------------------------------------- |  
    | LAT              | Latitude  (degrees)                               |
    | LON              | Longitude (degrees; from 0-360)                   |
    | TIME             | Timestamp (days since YYYY-N-D 0:0:0) (NetCDF Only) |
    | FH               | Forest canopy height (m)                          |
    | HREF             | Reference height above canopy (m)                 |
    | WS               | Wind speed at HREF (m/s), e.g., 10 m              |
    | CLU              | Forest clumping index (dimensionless)             |
    | LAI              | Leaf area index (m2/m2)                           |
    | VTYPE            | Vegetation type (dimensionless), VIIRS Only       |
    | FFRAC            | Forest fraction (dimensionless)                   |
    | UST              | Friction velocity (m/s)                           |
    | CSZ              | Cosine of the solar zenith angle (dimensionless)  |
    | Z0               | Total surface roughness length (m)                |
    | MOL              | Monin-Obukhov Length (m)                          |
    | FRP              | Total Fire Radiative Power (MW/grid cell area)    |

    **Table 2. Current User Namelist Options**

    | Namelist Option  | Namelist Description and Units                                                       |
    | ---------------  | ---------------------------------------------------------------------------------- |
    | infmt_opt        | integer for choosing 1D or 2D input file format (default = 0, 2D)                  |
    | nlat             | number of latitude cells (must match # of LAT in `file_vars` above)                |
    | nlon             | number of longitude cells (must match # of LON in `file_vars`above)                |
    | modlays          | number of model (below and above canopy) layers                                    |
    | modres           | above and below canopy model vertical resolution (m)                               |
    | ifcanwind        | logical canopy wind option (default = .FALSE.)                                     |
    | ifcanwaf         | logical canopy WAF option (default = .FALSE.)**                                    |
    | ifcaneddy        | logical canopy eddy Kz option (default = .FALSE.)                                  |
    | ifcanphot        | logical canopy photolysis option (default = .FALSE.)                               |
    | href_opt         | integer for using href_set in namelist (= 0) or array from file (= 1); default = 0 |
    | href_set         | user set real value of reference height above canopy associated with input wind speed (m) (only used if href_opt = 0)  |
    | z0ghc            | ratio of ground roughness length to canopy top height (Massman et al., 2017)       |
    | lamdars          | Value representing influence of roughness sublayer (Massman et al., 2017)          |
    | dx_opt           | 0=Calculation of dx resolution/distance from lon; 1=user set dx grid resolution    |
    | dx_set           | user set real value of grid resolution (m) only if dx_opt=1                        |
    | flameh_opt       | 0=Calculation of flame height from FRP (Byram, 1959); 1=user set flameh; 2=FRP calculation where available (active fires), otherwise user set flameh    |
    | flameh_set       | user set real value of flame height (m) only if flame_opt=1                        |
    | pai_opt          | integer (0=PAI fixed from Katul et al. 2004 veg types-->default; 1=PAI Massman et al. 2017 Eq. 19 calc; 2=PAI from model LAI+WAI; 3=user set PAI value) |
    | pai_set          | user set real value of PAI (default=4.0; only used if pai_opt=3)                   |
    | lu_opt           | integer (0=VIIRS LU type input mapped to Massman et al., currently only option )   |
    | z0_opt           | integer (0=use model input or 1 = vegtype dependent z0 for first estimate)         |
    | lai_thresh       | user set real value of LAI threshold for contiguous canopy (m2/m2)                 |
    | frt_thresh       | user set real value of forest fraction threshold for contiguous canopy             |
    | fch_thresh       | user set real value of canopy height  threshold for contiguous canopy (m)          |
    | rsl_opt          | user set option to include stability and Roughness SubLayer (RSL) effects  0 = off; 1 = on |

**Note:  If modres >> flameh then some error in WAF calculation will be incurred.  Suggestion is to use relative fine modres compared to average flame heights (at least <= 0.5 m) if WAF is required.


Main Citations (further references contained within):

Katul, G.G., Mahrt, L., Poggi, D., and Sanz, C. (2004). One- and two-equation models for canopy turbulence. Boundary-Layer Meteorol. 113: 81–109. https://doi.org/10.1023/B:BOUN.0000037333.48760.e5

Massman, W. J., J.M. Forthofer, and M.A. Finney. (2017). An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior. Canadian Journal of Forest Research. 47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354

Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone. Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243

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
