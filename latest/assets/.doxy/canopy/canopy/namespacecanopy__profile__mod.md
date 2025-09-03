

# Namespace canopy\_profile\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_profile\_mod**](namespacecanopy__profile__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_foliage**](#function-canopy_foliage) (MODLAYS MODLAYS, ZHC ZHC, ZCANMAX ZCANMAX, SIGMAU SIGMAU, SIGMA1 SIGMA1, FAFRACZINT FAFRACZINT) <br>_Compute canopy foliage and plant distribution functions._  |
|  subroutine | [**canopy\_zpd**](#function-canopy_zpd) (ZHC ZHC, FCLAI FCLAI, UBZREF UBZREF, Z0GHC Z0GHC, LAMBDARS LAMBDARS, CDRAG CDRAG, PAI PAI, FCH FCH, HREF HREF, Z0\_MOD Z0\_MOD, VTYPE VTYPE, LU\_OPT LU\_OPT, Z0\_OPT Z0\_OPT, d\_h d\_h, zo\_h zo\_h) <br>_Compute zero plane displacement height and surface roughness length._  |




























## Public Functions Documentation




### function canopy\_foliage 

_Compute canopy foliage and plant distribution functions._ 
```Fortran
subroutine canopy_profile_mod::canopy_foliage (
    MODLAYS MODLAYS,
    ZHC ZHC,
    ZCANMAX ZCANMAX,
    SIGMAU SIGMAU,
    SIGMA1 SIGMA1,
    FAFRACZINT FAFRACZINT
) 
```



Computes canopy foliage and plant distribution functions based on Massman et al. (2017) algorithms. Calculates foliage area fraction distribution throughout the canopy layers. 

**Parameters:**


* `MODLAYS` Number of model layers 
* `ZHC` Height of each canopy model layer (z/h) 
* `ZCANMAX` Height of maximum foliage area density (z/h) 
* `SIGMAU` Standard deviation of shape function above zcanmax (z/h) 
* `SIGMA1` Standard deviation of shape function below zcanmax (z/h) 
* `FAFRACZINT` Foliage area fraction at each canopy model layer [output] 



**Author:**

P.C. Campbell 




**Date:**

October 2022 




**Note:**

Based on Massman et al. (2017) algorithms 





        

<hr>



### function canopy\_zpd 

_Compute zero plane displacement height and surface roughness length._ 
```Fortran
subroutine canopy_profile_mod::canopy_zpd (
    ZHC ZHC,
    FCLAI FCLAI,
    UBZREF UBZREF,
    Z0GHC Z0GHC,
    LAMBDARS LAMBDARS,
    CDRAG CDRAG,
    PAI PAI,
    FCH FCH,
    HREF HREF,
    Z0_MOD Z0_MOD,
    VTYPE VTYPE,
    LU_OPT LU_OPT,
    Z0_OPT Z0_OPT,
    d_h d_h,
    zo_h zo_h
) 
```



Computes zero plane displacement height and total surface roughness length based on Massman et al. (2017) algorithms. Uses foliage characteristics, drag coefficient, and in-canopy plant distribution functions. 

**Parameters:**


* `ZHC` Height of each canopy model layer (z/h) 
* `FCLAI` Foliage area fraction at each canopy model layer 
* `UBZREF` Mean wind speed at reference height (m/s) 
* `Z0GHC` Ground roughness length to canopy height ratio 
* `LAMBDARS` Leaf-scale resistance adjustment factor 
* `CDRAG` Drag coefficient (nondimensional) 
* `PAI` Plant/foliage area index (nondimensional) 
* `FCH` Grid cell canopy height (m) 
* `HREF` Reference height (m) 
* `Z0_MOD` Surface roughness length modifier 
* `VTYPE` Grid cell dominant vegetation type 
* `LU_OPT` Integer for LU type from model 
* `Z0_OPT` Integer for roughness length option 
* `d_h` Zero plane displacement height normalized by canopy height [output] 
* `zo_h` Surface roughness length normalized by canopy height [output] 



**Author:**

P.C. Campbell 




**Date:**

June 2022 




**Note:**

Based on Massman et al. (2017) algorithms 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_profile_mod.F90`

