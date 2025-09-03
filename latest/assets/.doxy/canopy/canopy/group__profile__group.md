

# Group profile\_group



[**Modules**](modules.md) **>** [**profile\_group**](group__profile__group.md)



_Routines for calculating canopy profiles and foliage distributions._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**cdrag**](#variable-cdrag)  <br>_Drag coefficient (nondimensional)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fch**](#variable-fch)  <br>_Grid cell canopy height (m)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ffrac**](#variable-ffrac)  <br>_Grid cell forest fraction._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**firetype**](#variable-firetype)  <br>_1 = Above Canopy Fire; 0 = Below Canopy Fire; -1 No Canopy_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lai**](#variable-lai)  <br>_Grid cell leaf area index._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lu\_opt**](#variable-lu_opt)  <br>_Integer for LU type from model mapped to Massman et al._  |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**pai**](#variable-pai)  <br>_Plant/foliage area index (nondimensional)_  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pai\_opt**](#variable-pai_opt)  <br>_Integer for PAI values used or calculated (default = 0)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pai\_set**](#variable-pai_set)  <br>_Real value for PAI set values used (default = 4.0)_  |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**sigma1**](#variable-sigma1)  <br>_Standard deviation of shape function below zcanmax (z/h)_  |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**sigmau**](#variable-sigmau)  <br>_Standard deviation of shape function above zcanmax (z/h)_  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vtype**](#variable-vtype)  <br>_Grid cell dominant vegetation type._  |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**zcanmax**](#variable-zcanmax)  <br>_Height of maximum foliage area density (z/h) (nondimensional)_  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_profile\_mod::canopy\_foliage**](#function-canopy_foliage) (MODLAYS MODLAYS, ZHC ZHC, ZCANMAX ZCANMAX, SIGMAU SIGMAU, SIGMA1 SIGMA1, FAFRACZINT FAFRACZINT) <br>_Compute canopy foliage and plant distribution functions._  |




























## Detailed Description


This group contains subroutines for calculating canopy parameters, foliage area density distributions, and zero-plane displacement based on Massman et al. (2017) algorithms and vegetation characteristics. 


    
## Public Attributes Documentation




### variable cdrag 

_Drag coefficient (nondimensional)_ 
```Fortran
real(rk), intent(out) cdrag;
```




<hr>



### variable fch 

_Grid cell canopy height (m)_ 
```Fortran
real(rk), intent(in) fch;
```




<hr>



### variable ffrac 

_Grid cell forest fraction._ 
```Fortran
real(rk), intent(in) ffrac;
```




<hr>



### variable firetype 

_1 = Above Canopy Fire; 0 = Below Canopy Fire; -1 No Canopy_ 
```Fortran
integer, intent(out) firetype;
```




<hr>



### variable lai 

_Grid cell leaf area index._ 
```Fortran
real(rk), intent(in) lai;
```




<hr>



### variable lu\_opt 

_Integer for LU type from model mapped to Massman et al._ 
```Fortran
integer, intent(in) lu_opt;
```



(default = 0/VIIRS) 


        

<hr>



### variable pai 

_Plant/foliage area index (nondimensional)_ 
```Fortran
real(rk), intent(out) pai;
```




<hr>



### variable pai\_opt 

_Integer for PAI values used or calculated (default = 0)_ 
```Fortran
integer, intent(in) pai_opt;
```




<hr>



### variable pai\_set 

_Real value for PAI set values used (default = 4.0)_ 
```Fortran
real(rk), intent(in) pai_set;
```




<hr>



### variable sigma1 

_Standard deviation of shape function below zcanmax (z/h)_ 
```Fortran
real(rk), intent(out) sigma1;
```




<hr>



### variable sigmau 

_Standard deviation of shape function above zcanmax (z/h)_ 
```Fortran
real(rk), intent(out) sigmau;
```




<hr>



### variable vtype 

_Grid cell dominant vegetation type._ 
```Fortran
integer, intent(in) vtype;
```




<hr>



### variable zcanmax 

_Height of maximum foliage area density (z/h) (nondimensional)_ 
```Fortran
real(rk), intent(out) zcanmax;
```




<hr>
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

------------------------------


