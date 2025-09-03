

# Group var3din\_inputs



[**Modules**](modules.md) **>** [**var3din\_inputs**](group__var3din__inputs.md)



_Compute integral of incremental fractional foliage shape function from GEDI PAVD._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**above**](#variable-above)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**area**](#variable-area)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**associated**](#variable-associated)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**below**](#variable-below)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**canopy**](#variable-canopy)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cell**](#variable-cell)  <br> |
|  real(rk), dimension(nondimensional), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**density**](#variable-density)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**deviation**](#variable-deviation)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fch**](#variable-fch)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**foliage**](#variable-foliage)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**from**](#variable-from)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**function**](#variable-function)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**gedi**](#variable-gedi)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**grid**](#variable-grid)  <br> |
|  real(rk), dimension(dimensionless), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**h**](#variable-h)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**height**](#variable-height)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**heights**](#variable-heights)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**level**](#variable-level)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**maximum**](#variable-maximum)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**mid**](#variable-mid)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**of**](#variable-of)  <br> |
|  real(rk), dimension([**m**](canopy__bioemi__mod_8F90.md#variable-m)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pavd**](#variable-pavd)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pavd\_in**](#variable-pavd_in)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pavd\_levs**](#variable-pavd_levs)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**plant**](#variable-plant)  <br> |
|  real(rk), dimension(m²/m³), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**profile**](#variable-profile)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)), allocatable | [**shape**](#variable-shape)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**sigma1**](#variable-sigma1)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**sigmau**](#variable-sigmau)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**standard**](#variable-standard)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**volume**](#variable-volume)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**with**](#variable-with)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**z**](#variable-z)  <br> |
|  real(rk), dimension([**z**](canopy__var3din__mod_8F90.md#variable-z)/[**h**](canopy__var3din__mod_8F90.md#variable-h)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zcanmax**](#variable-zcanmax)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zcanmax\_in**](#variable-zcanmax_in)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zhc**](#variable-zhc)  <br> |












































## Detailed Description


This subroutine converts interpolated GEDI 3D PAVD (Plant Area Volume Density) profiles into fractional foliage shape functions using the algorithms from Massman et al. (2017). The process includes:
* Interpolating input PAVD data to the canopy model vertical resolution
* Determining the height of maximum foliage area density (ZCANMAX) from observed PAVD
* Calculating incremental foliage shape functions using Gaussian distributions
* Computing fractional cumulative foliage distributions
* Integrating the foliage shape functions for canopy structure parameterization






**Parameters:**


* `ZCANMAX_IN` Input height of maximum foliage area density (z/h) (nondimensional) 
* `SIGMAU` Standard deviation of shape function above zcanmax (z/h) 
* `SIGMA1` Standard deviation of shape function below zcanmax (z/h) 
* `FCH` Grid cell canopy height (m) from GEDI 
* `ZHC` Dimensionless height coordinate (z/h) 
* `PAVD_IN` Plant Area Volume Density profile (m²/m³) 
* `PAVD_LEVS` Associated mid-level heights for PAVD data (m) 
* `FAFRACZINT` Integral of incremental fractional foliage shape function 




    
## Public Attributes Documentation




### variable above 

```Fortran
real(rk), intent(in) above;
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



### variable input 

```Fortran
real(rk), intent(in) input;
```




<hr>



### variable level 

```Fortran
real(rk), intent(in) level;
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

------------------------------


