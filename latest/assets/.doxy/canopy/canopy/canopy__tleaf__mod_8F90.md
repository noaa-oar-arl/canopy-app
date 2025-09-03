

# File canopy\_tleaf\_mod.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_tleaf\_mod.F90**](canopy__tleaf__mod_8F90.md)

[Go to the source code of this file](canopy__tleaf__mod_8F90_source.md)

_Leaf temperature calculation module for canopy model._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**canopy\_tleaf\_mod**](namespacecanopy__tleaf__mod.md) <br> |








## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), parameter | [**atemp\_1\_shade**](#variable-atemp_1_shade)   = `-12.846\_rk`<br> |
|  real(rk), parameter | [**atemp\_1\_sun**](#variable-atemp_1_sun)   = `-13.891\_rk`<br> |
|  real(rk), parameter | [**atemp\_2\_shade**](#variable-atemp_2_shade)   = `-11.343\_rk`<br> |
|  real(rk), parameter | [**atemp\_2\_sun**](#variable-atemp_2_sun)   = `-12.322\_rk`<br> |
|  real(rk), parameter | [**atemp\_3\_shade**](#variable-atemp_3_shade)   = `-1.068\_rk`<br> |
|  real(rk), parameter | [**atemp\_3\_sun**](#variable-atemp_3_sun)   = `-1.032\_rk`<br> |
|  real(rk), parameter | [**atemp\_4\_shade**](#variable-atemp_4_shade)   = `-5.551\_rk`<br> |
|  real(rk), parameter | [**atemp\_4\_sun**](#variable-atemp_4_sun)   = `-5.172\_rk`<br> |
|  real(rk), parameter | [**atemp\_5\_shade**](#variable-atemp_5_shade)   = `-5.955\_rk`<br> |
|  real(rk), parameter | [**atemp\_5\_sun**](#variable-atemp_5_sun)   = `-5.589\_rk`<br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**atemp\_shade**](#variable-atemp_shade)  <br>_Regression coefficient A for shade leaves._  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**atemp\_sun**](#variable-atemp_sun)  <br>_Regression coefficient A for sun leaves (Silva et al., 2020)_  |
|  real(rk), parameter | [**btemp\_1\_shade**](#variable-btemp_1_shade)   = `1.060\_rk`<br> |
|  real(rk), parameter | [**btemp\_1\_sun**](#variable-btemp_1_sun)   = `1.064\_rk`<br> |
|  real(rk), parameter | [**btemp\_2\_shade**](#variable-btemp_2_shade)   = `1.053\_rk`<br> |
|  real(rk), parameter | [**btemp\_2\_sun**](#variable-btemp_2_sun)   = `1.057\_rk`<br> |
|  real(rk), parameter | [**btemp\_3\_shade**](#variable-btemp_3_shade)   = `1.031\_rk`<br> |
|  real(rk), parameter | [**btemp\_3\_sun**](#variable-btemp_3_sun)   = `1.031\_rk`<br> |
|  real(rk), parameter | [**btemp\_4\_shade**](#variable-btemp_4_shade)   = `1.051\_rk`<br> |
|  real(rk), parameter | [**btemp\_4\_sun**](#variable-btemp_4_sun)   = `1.050\_rk`<br> |
|  real(rk), parameter | [**btemp\_5\_shade**](#variable-btemp_5_shade)   = `1.053\_rk`<br> |
|  real(rk), parameter | [**btemp\_5\_sun**](#variable-btemp_5_sun)   = `1.051\_rk`<br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**btemp\_shade**](#variable-btemp_shade)  <br>_Regression coefficient B for shade leaves._  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**btemp\_sun**](#variable-btemp_sun)  <br>_Regression coefficient B for sun leaves._  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fch**](#variable-fch)  <br>_Model input canopy height (m)_  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fsun**](#variable-fsun)  <br>_Sunlit/Shaded fraction from photolysis correction factor._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**i**](#variable-i)  <br>_Loop index._  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**temp2**](#variable-temp2)  <br>_Model input 2-m Temperature (K)_  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**tleaf\_ave**](#variable-tleaf_ave)  <br>_Ave Leaf temp for sun/shaded leaves (K)_  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**tleaf\_shade**](#variable-tleaf_shade)  <br>_Leaf temp for shaded leaves (K)_  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**tleaf\_sun**](#variable-tleaf_sun)  <br>_Leaf temp for sunlit leaves (K)_  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zk**](#variable-zk)  <br>_Input model heights (m)_  |












































## Detailed Description


This module implements linear interpolation methods for calculating leaf temperatures for sunlit and shaded leaves throughout the canopy based on Silva et al. (2020) algorithms. 

**Author:**

P.C. Campbell 




**Date:**

June 2023 




**Version:**

1.0 





    
## Public Attributes Documentation




### variable atemp\_1\_shade 

```Fortran
real(rk), parameter atemp_1_shade;
```




<hr>



### variable atemp\_1\_sun 

```Fortran
real(rk), parameter atemp_1_sun;
```




<hr>



### variable atemp\_2\_shade 

```Fortran
real(rk), parameter atemp_2_shade;
```




<hr>



### variable atemp\_2\_sun 

```Fortran
real(rk), parameter atemp_2_sun;
```




<hr>



### variable atemp\_3\_shade 

```Fortran
real(rk), parameter atemp_3_shade;
```




<hr>



### variable atemp\_3\_sun 

```Fortran
real(rk), parameter atemp_3_sun;
```




<hr>



### variable atemp\_4\_shade 

```Fortran
real(rk), parameter atemp_4_shade;
```




<hr>



### variable atemp\_4\_sun 

```Fortran
real(rk), parameter atemp_4_sun;
```




<hr>



### variable atemp\_5\_shade 

```Fortran
real(rk), parameter atemp_5_shade;
```




<hr>



### variable atemp\_5\_sun 

```Fortran
real(rk), parameter atemp_5_sun;
```




<hr>



### variable atemp\_shade 

_Regression coefficient A for shade leaves._ 
```Fortran
real(rk), dimension(size(zk)) atemp_shade;
```




<hr>



### variable atemp\_sun 

_Regression coefficient A for sun leaves (Silva et al., 2020)_ 
```Fortran
real(rk), dimension(size(zk)) atemp_sun;
```




<hr>



### variable btemp\_1\_shade 

```Fortran
real(rk), parameter btemp_1_shade;
```




<hr>



### variable btemp\_1\_sun 

```Fortran
real(rk), parameter btemp_1_sun;
```




<hr>



### variable btemp\_2\_shade 

```Fortran
real(rk), parameter btemp_2_shade;
```




<hr>



### variable btemp\_2\_sun 

```Fortran
real(rk), parameter btemp_2_sun;
```




<hr>



### variable btemp\_3\_shade 

```Fortran
real(rk), parameter btemp_3_shade;
```




<hr>



### variable btemp\_3\_sun 

```Fortran
real(rk), parameter btemp_3_sun;
```




<hr>



### variable btemp\_4\_shade 

```Fortran
real(rk), parameter btemp_4_shade;
```




<hr>



### variable btemp\_4\_sun 

```Fortran
real(rk), parameter btemp_4_sun;
```




<hr>



### variable btemp\_5\_shade 

```Fortran
real(rk), parameter btemp_5_shade;
```




<hr>



### variable btemp\_5\_sun 

```Fortran
real(rk), parameter btemp_5_sun;
```




<hr>



### variable btemp\_shade 

_Regression coefficient B for shade leaves._ 
```Fortran
real(rk), dimension(size(zk)) btemp_shade;
```




<hr>



### variable btemp\_sun 

_Regression coefficient B for sun leaves._ 
```Fortran
real(rk), dimension(size(zk)) btemp_sun;
```




<hr>



### variable fch 

_Model input canopy height (m)_ 
```Fortran
real(rk), intent(in) fch;
```




<hr>



### variable fsun 

_Sunlit/Shaded fraction from photolysis correction factor._ 
```Fortran
real(rk), dimension(:), intent(in) fsun;
```




<hr>



### variable i 

_Loop index._ 
```Fortran
integer i;
```



Calculate linear change in parameters interpolated to Silva et al. 5 layer canopy regions Above canopy, Tleaf = Tair Level 1 - 2 Level 2 - 3 Level 3 - 4 Level 4 - Bottom 


        

<hr>



### variable temp2 

_Model input 2-m Temperature (K)_ 
```Fortran
real(rk), intent(in) temp2;
```




<hr>



### variable tleaf\_ave 

_Ave Leaf temp for sun/shaded leaves (K)_ 
```Fortran
real(rk), dimension(size(zk)), intent(out) tleaf_ave;
```




<hr>



### variable tleaf\_shade 

_Leaf temp for shaded leaves (K)_ 
```Fortran
real(rk), dimension(size(zk)), intent(out) tleaf_shade;
```




<hr>



### variable tleaf\_sun 

_Leaf temp for sunlit leaves (K)_ 
```Fortran
real(rk), dimension(size(zk)), intent(out) tleaf_sun;
```




<hr>



### variable zk 

_Input model heights (m)_ 
```Fortran
real(rk), dimension(:), intent(in) zk;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_tleaf_mod.F90`

