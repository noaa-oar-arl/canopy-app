

# Group canopy\_wind



[**Modules**](modules.md) **>** [**canopy\_wind**](group__canopy__wind.md)



_Wind speed profile calculations within and above canopy._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_wind\_mod::canopy\_wind\_most**](#function-canopy_wind_most) (HCM HCM, ZK ZK, FAFRACK FAFRACK, UBZREF UBZREF, Z0GHC Z0GHC, CDRAG CDRAG, PAI PAI, HREF HREF, D\_H D\_H, ZO\_H ZO\_H, LAMBDARS LAMBDARS, CANBOT\_OUT CANBOT\_OUT, CANTOP\_OUT CANTOP\_OUT, CANWIND CANWIND) <br>_Calculate canopy wind speed profile using MOST._  |




























## Public Functions Documentation




### function canopy\_wind\_most 

_Calculate canopy wind speed profile using MOST._ 
```Fortran
subroutine canopy_wind_mod::canopy_wind_most (
    HCM HCM,
    ZK ZK,
    FAFRACK FAFRACK,
    UBZREF UBZREF,
    Z0GHC Z0GHC,
    CDRAG CDRAG,
    PAI PAI,
    HREF HREF,
    D_H D_H,
    ZO_H ZO_H,
    LAMBDARS LAMBDARS,
    CANBOT_OUT CANBOT_OUT,
    CANTOP_OUT CANTOP_OUT,
    CANWIND CANWIND
) 
```



Computes mean wind speed for given height (z) below the canopy top using Monin-Obukhov Similarity Theory above canopy and Massman et al. (2017) parameterization within canopy




**Parameters:**


* `HCM` Height of canopy top (m) 
* `ZK` Above/Below canopy height, z (m) 
* `FAFRACK` Fractional (z) shapes of the plant surface distribution (dimensionless) 
* `UBZREF` Mean wind speed at reference height (ZREF) (m/s) 
* `Z0GHC` Ratio of ground roughness length to canopy top height (dimensionless) 
* `CDRAG` Drag coefficient (dimensionless) 
* `PAI` Total plant/foliage area index (dimensionless) 
* `HREF` Reference Height above the canopy (ZREF = HCM + HREF; m) 
* `D_H` Zero plane displacement heights (dimensionless) 
* `ZO_H` Surface (soil+veg) roughness lengths (dimensionless) 
* `LAMBDARS` Value representing influence of roughness sublayer (dimensionless) 
* `CANBOT_OUT` Canopy bottom wind reduction factor = canbot (dimensionless) 
* `CANTOP_OUT` Canopy top wind reduction factor = cantop (dimensionless) 
* `CANWIND` Mean canopy wind speed at current z (m/s) 



**Author:**

P. C. Campbell 




**Date:**

Jun 2022 




**Note:**

Based on Massman et al. (2017) algorithms below the canopy and use of MOST above canopy Massman, W.J., J.M. Forthofer, M.A. Finney (2017). An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior. Canadian Journal of Forest Research, 47(5), 594-599. [https://doi.org/10.1139/cjfr-2016-0354](https://doi.org/10.1139/cjfr-2016-0354) 





        

<hr>

------------------------------


