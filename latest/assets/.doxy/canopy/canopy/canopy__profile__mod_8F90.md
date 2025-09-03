

# File canopy\_profile\_mod.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_profile\_mod.F90**](canopy__profile__mod_8F90.md)

[Go to the source code of this file](canopy__profile__mod_8F90_source.md)

_Canopy profile and foliage distribution calculations._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**canopy\_profile\_mod**](namespacecanopy__profile__mod.md) <br> |








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












































## Detailed Description


This module implements canopy profile calculations including canopy parameters, foliage distribution, and zero-plane displacement calculations based on vegetation characteristics. 

**Author:**

P.C. Campbell 




**Date:**

June 2022 




**Version:**

1.0 





    
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

------------------------------
The documentation for this class was generated from the following file `src/canopy_profile_mod.F90`

