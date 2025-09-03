

# File canopy\_var3din\_mod.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_var3din\_mod.F90**](canopy__var3din__mod_8F90.md)

[Go to the source code of this file](canopy__var3din__mod_8F90_source.md)

_3D Variable Input Module_ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**canopy\_var3din\_mod**](namespacecanopy__var3din__mod.md) <br> |








## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**above**](#variable-above)  <br> |
|  real(rk) | [**actual**](#variable-actual)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**area**](#variable-area)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**associated**](#variable-associated)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**below**](#variable-below)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**canopy**](#variable-canopy)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cell**](#variable-cell)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)) | [**coordinates**](#variable-coordinates)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**counters**](#variable-counters)  <br> |
|  real(rk), dimension(nondimensional), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**density**](#variable-density)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**deviation**](#variable-deviation)  <br> |
|  real(rk), dimension(:), allocatable | [**fafracz**](#variable-fafracz)  <br> |
|  real(rk), dimension(:), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**fafraczint**](#variable-fafraczint)  <br> |
|  real(rk), dimension(:), allocatable | [**fainc**](#variable-fainc)  <br> |
|  real(rk) | [**fatot**](#variable-fatot)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fch**](#variable-fch)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**foliage**](#variable-foliage)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**fractional**](#variable-fractional)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**from**](#variable-from)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**function**](#variable-function)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**gedi**](#variable-gedi)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**grid**](#variable-grid)  <br> |
|  real(rk), dimension(dimensionless), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**h**](#variable-h)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**height**](#variable-height)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**heights**](#variable-heights)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**i**](#variable-i)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**incremental**](#variable-incremental)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**integral**](#variable-integral)  <br> |
|  real(rk) | [**interpolated**](#variable-interpolated)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**lev**](#variable-lev)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**level**](#variable-level)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**loop**](#variable-loop)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**maximum**](#variable-maximum)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**mid**](#variable-mid)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**of**](#variable-of)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pavd**](#variable-pavd)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pavd\_in**](#variable-pavd_in)  <br> |
|  real(rk), dimension(size([**zhc**](canopy__var3din__mod_8F90.md#variable-zhc))) | [**pavd\_interp**](#variable-pavd_interp)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pavd\_levs**](#variable-pavd_levs)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**plant**](#variable-plant)  <br> |
|  real(rk), dimension(m²/m³), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**profile**](#variable-profile)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**shape**](#variable-shape)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**sigma1**](#variable-sigma1)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**sigmau**](#variable-sigmau)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**standard**](#variable-standard)  <br> |
|  real(rk) | [**total**](#variable-total)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**volume**](#variable-volume)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**with**](#variable-with)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**z**](#variable-z)  <br> |
|  real(rk), dimension([**z**](canopy__var3din__mod_8F90.md#variable-z)/[**h**](canopy__var3din__mod_8F90.md#variable-h)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zcanmax**](#variable-zcanmax)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zcanmax\_in**](#variable-zcanmax_in)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zhc**](#variable-zhc)  <br> |
|  real(rk), dimension(size([**zhc**](canopy__var3din__mod_8F90.md#variable-zhc))) | [**zk**](#variable-zk)  <br> |












































## Detailed Description


This module contains subroutines for processing 3D variable inputs, particularly for converting GEDI PAVD (Plant Area Volume Density) profiles into fractional foliage shape functions. The module handles interpolation of observed PAVD profiles to user-defined canopy model resolutions.




**Author:**

Patrick C. Campbell 




**Date:**

July 2023


\references Massman, W.J., Forthofer, J.M., and Finney, M.A.: An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior. Canadian Journal of Forest Research. 47(5): 594-603. [https://doi.org/10.1139/cjfr-2016-0354](https://doi.org/10.1139/cjfr-2016-0354) 


    
## Public Attributes Documentation




### variable above 

```Fortran
real(rk), intent(in) above;
```




<hr>



### variable actual 

```Fortran
real(rk) actual;
```




<hr>



### variable area 

```Fortran
real(rk), intent(in) area;
```




<hr>



### variable associated 

```Fortran
real(rk), intent(in) associated;
```




<hr>



### variable below 

```Fortran
real(rk), intent(in) below;
```




<hr>



### variable canopy 

```Fortran
real(rk), intent(in) canopy;
```




<hr>



### variable cell 

```Fortran
real(rk), intent(in) cell;
```




<hr>



### variable coordinates 

```Fortran
real(rk), dimension (m) coordinates;
```




<hr>



### variable counters 

```Fortran
integer counters;
```




<hr>



### variable density 

```Fortran
real(rk), dimension (nondimensional), intent(in) density;
```




<hr>



### variable deviation 

```Fortran
real(rk), intent(in) deviation;
```




<hr>



### variable fafracz 

```Fortran
real(rk), dimension(:), allocatable fafracz;
```




<hr>



### variable fafraczint 

```Fortran
real(rk), dimension(:), intent(out) fafraczint;
```




<hr>



### variable fainc 

```Fortran
real(rk), dimension(:), allocatable fainc;
```




<hr>



### variable fatot 

```Fortran
real(rk) fatot;
```




<hr>



### variable fch 

```Fortran
real(rk), intent(in) fch;
```




<hr>



### variable foliage 

```Fortran
real(rk), intent(out), allocatable foliage;
```




<hr>



### variable fractional 

```Fortran
real(rk), intent(out), allocatable fractional;
```




<hr>



### variable from 

```Fortran
real(rk), intent(in) from;
```




<hr>



### variable function 

```Fortran
real(rk), intent(out), allocatable function;
```




<hr>



### variable gedi 

```Fortran
real(rk), intent(in) gedi;
```




<hr>



### variable grid 

```Fortran
real(rk), intent(in) grid;
```




<hr>



### variable h 

```Fortran
real(rk), dimension (dimensionless), intent(in) h;
```




<hr>



### variable height 

```Fortran
real(rk), dimension (m), intent(in) height;
```




<hr>



### variable heights 

```Fortran
real(rk), intent(in) heights;
```




<hr>



### variable i 

```Fortran
integer i;
```




<hr>



### variable incremental 

```Fortran
real(rk), intent(out), allocatable incremental;
```




<hr>



### variable input 

```Fortran
real(rk), intent(in) input;
```




<hr>



### variable integral 

```Fortran
real(rk), intent(out) integral;
```




<hr>



### variable interpolated 

```Fortran
real(rk) interpolated;
```




<hr>



### variable lev 

```Fortran
integer lev;
```




<hr>



### variable level 

```Fortran
real(rk), intent(in) level;
```




<hr>



### variable loop 

```Fortran
integer loop;
```




<hr>



### variable maximum 

```Fortran
real(rk), intent(in) maximum;
```




<hr>



### variable mid 

```Fortran
real(rk), intent(in) mid;
```




<hr>



### variable of 

```Fortran
real(rk), intent(out) of;
```




<hr>



### variable pavd 

```Fortran
real(rk), dimension (m), intent(in) pavd;
```




<hr>



### variable pavd\_in 

```Fortran
real(rk), dimension(:), intent(in) pavd_in;
```




<hr>



### variable pavd\_interp 

```Fortran
real(rk), dimension(size(zhc)) pavd_interp;
```




<hr>



### variable pavd\_levs 

```Fortran
real(rk), dimension(:), intent(in) pavd_levs;
```




<hr>



### variable plant 

```Fortran
real(rk), intent(in) plant;
```




<hr>



### variable profile 

```Fortran
real(rk), dimension (m²/m³), intent(in) profile;
```




<hr>



### variable shape 

```Fortran
real(rk), intent(out), allocatable shape;
```




<hr>



### variable sigma1 

```Fortran
real(rk), intent(in) sigma1;
```




<hr>



### variable sigmau 

```Fortran
real(rk), intent(in) sigmau;
```




<hr>



### variable standard 

```Fortran
real(rk), intent(in) standard;
```




<hr>



### variable total 

```Fortran
real(rk) total;
```




<hr>



### variable volume 

```Fortran
real(rk), intent(in) volume;
```




<hr>



### variable with 

```Fortran
real(rk), intent(in) with;
```




<hr>



### variable z 

```Fortran
real(rk), intent(in) z;
```




<hr>



### variable zcanmax 

```Fortran
real(rk), dimension (z/h), intent(in) zcanmax;
```




<hr>



### variable zcanmax\_in 

```Fortran
real(rk), intent(in) zcanmax_in;
```




<hr>



### variable zhc 

```Fortran
real(rk), dimension(:), intent(in) zhc;
```




<hr>



### variable zk 

```Fortran
real(rk), dimension(size(zhc)) zk;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_var3din_mod.F90`

