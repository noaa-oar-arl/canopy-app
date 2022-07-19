# canopy-app
Repository for low-level, stand-alone/column canopy parameterizations for testing and application to gridded atmospheric composition/air quality models.

Author(s):

Patrick Campbell, Zachary Moon, and Wei-Ting Hung

Compile, edit namelist, and run canopy model:
- `make`
- `namelist.canopy`
- `./canopy`

Canopy is parameterized by foliage distribution shape functions and parameters for different vegetation types (Katul et al., 2004; Massman et al., 2017)

- `canopy_parm_mod.F90`

Current Canopy-App components:

1.  In-Canopy Winds and Wind Adjustment Factor (WAF) for wildland fire spread model applications (Massman et al., 2017).

- `canopy_wind_mod.F90`
- `canopy_waf_mod.F90`

2.  In-Canopy radiative transfer (i.e., photolysis rate adjusment) for air quality applications (in progress...)  (Makar et al., 2017)

- `canopy_prad_mod.F90`

3.  In-Canopy vertical diffusion (i.e., eddy diffusivity adjustment) for air quality applications (in progress...) (Makar et al., 2017)

- `canopy_vdiff_mod.F90`


Citations:

Katul, G.G., Mahrt, L., Poggi, D., and Sanz, C. (2004). One- and two-equation models for canopy turbulence. Boundary-Layer Meteorol. 113: 81–109. doi:10.1023/B:BOUN.0000037333.48760.e5.

Massman, W. J., J.M. Forthofer, and M.A. Finney. (2017). An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior. Canadian Journal of Forest Research. 47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354

Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone. Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243


## Development

To set up the pre-commit hooks, invoke
```
pre-commit install --install-hooks
```
after [installing `pre-commit`](https://pre-commit.com/#installation).
