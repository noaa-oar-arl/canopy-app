# canopy_wind-fire
Stand-alone/column canopy winds and Wind Adjustment Factor (WAF) for wildland fire spread model applications.  

Model based on Massman et al. (2017).

Author:

Patrick Campbell

Usage:

Adjust canopy/vegetation type parameters starting on Line 30 of canopy_driver.F90 for desired example (or input based on wrapper)

Compile and run canopy model:
- `make`
- `./canopy`

Citation:

W.J. Massman, J.M. Forthofer, and M.A. Finney. An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior. Canadian Journal of Forest Research. 47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354
