

# Group canopy\_fire



[**Modules**](modules.md) **>** [**canopy\_fire**](group__canopy__fire.md)



_Fire-related calculations for flame height and wind adjustment factors._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_fire\_mod::canopy\_flameh**](#function-canopy_flameh) (FLAMEH\_OPT FLAMEH\_OPT, FLAMEH\_SET FLAMEH\_SET, DX DX, MODRES MODRES, FRP\_IN FRP\_IN, FRP\_FAC FRP\_FAC, FCH FCH, LU\_OPT LU\_OPT, VTYPE VTYPE, FLAMEH\_CAL FLAMEH\_CAL, MIDFLAMEPOINT MIDFLAMEPOINT, FLAMEH FLAMEH) <br>_Calculate flame height and midflame point._  |
|  subroutine | [**canopy\_fire\_mod::canopy\_waf**](#function-canopy_waf) (HCM HCM, LAMBDARS LAMBDARS, HREF HREF, FLAMEH FLAMEH, FIRETYPE FIRETYPE, D\_H D\_H, ZO\_H ZO\_H, CANBOTMID CANBOTMID, CANTOPMID CANTOPMID, WAF WAF) <br>_Calculate Wind Adjustment Factor for fire spread._  |




























## Public Functions Documentation




### function canopy\_flameh 

_Calculate flame height and midflame point._ 
```Fortran
subroutine canopy_fire_mod::canopy_flameh (
    FLAMEH_OPT FLAMEH_OPT,
    FLAMEH_SET FLAMEH_SET,
    DX DX,
    MODRES MODRES,
    FRP_IN FRP_IN,
    FRP_FAC FRP_FAC,
    FCH FCH,
    LU_OPT LU_OPT,
    VTYPE VTYPE,
    FLAMEH_CAL FLAMEH_CAL,
    MIDFLAMEPOINT MIDFLAMEPOINT,
    FLAMEH FLAMEH
) 
```



Computes flame height and midflamepoint layer needed for WAF calculation using various methods including FRP-based calculations and user-defined values




**Parameters:**


* `FLAMEH_OPT` Integer for flame height calculation option (0=calculate, 1=user set, etc.) 
* `FLAMEH_SET` User set flame height value (m) 
* `DX` Grid cell distance using haversine formula (m) 
* `MODRES` Canopy model input vertical resolution (m) 
* `FRP_IN` Model input Fire Radiative Power (MW/grid cell area) 
* `FRP_FAC` FRP tuning factor for flame height calculation 
* `FCH` Grid cell canopy height (m) 
* `LU_OPT` Supported land use classifications 
* `VTYPE` Dominant vegetation type 
* `FLAMEH_CAL` Option of vegetation type dependent FRP to flame height relationships 
* `MIDFLAMEPOINT` Index of the mid-flame point 
* `FLAMEH` Calculated flame height (m) 



**Author:**

P. C. Campbell 




**Date:**

Oct 2022 




**Note:**

Supports multiple flame height calculation methods and crown fire detection 





        

<hr>



### function canopy\_waf 

_Calculate Wind Adjustment Factor for fire spread._ 
```Fortran
subroutine canopy_fire_mod::canopy_waf (
    HCM HCM,
    LAMBDARS LAMBDARS,
    HREF HREF,
    FLAMEH FLAMEH,
    FIRETYPE FIRETYPE,
    D_H D_H,
    ZO_H ZO_H,
    CANBOTMID CANBOTMID,
    CANTOPMID CANTOPMID,
    WAF WAF
) 
```



Computes Wind Adjustment Factor for fire spread for either sub-canopy or above-canopy fires using Massman et al. (2017) formulations




**Parameters:**


* `HCM` Height of canopy top (m) 
* `LAMBDARS` Value representing influence of roughness sublayer (dimensionless) 
* `HREF` Reference height above the canopy (m) 
* `FLAMEH` Flame height for above canopy fire (m) 
* `FIRETYPE` Fire type: 1=Above Canopy Fire, 0=Below Canopy Fire 
* `CANBOTMID` Mid-flame canopy bottom wind reduction factor (dimensionless) 
* `CANTOPMID` Mid-flame canopy top wind reduction factor (dimensionless) 
* `D_H` Zero-plane displacement height (d/h) (dimensionless) 
* `ZO_H` Surface (soil+veg) roughness length (zo/h) (dimensionless) 
* `WAF` Wind Adjustment Factor (dimensionless) 



**Author:**

P. C. Campbell 




**Date:**

Jun 2022 




**Note:**

Based on Massman et al. (2017) algorithms for fire behavior prediction Massman, W.J., J.M. Forthofer, M.A. Finney (2017). An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior. Canadian Journal of Forest Research, 47(5), 594-599. [https://doi.org/10.1139/cjfr-2016-0354](https://doi.org/10.1139/cjfr-2016-0354) 





        

<hr>

------------------------------


