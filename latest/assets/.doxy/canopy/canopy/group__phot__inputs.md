

# Group phot\_inputs



[**Modules**](modules.md) **>** [**phot\_inputs**](group__phot__inputs.md)



_Calculate photolysis attenuation in canopy._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**angle**](#variable-angle)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**area**](#variable-area)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**clu**](#variable-clu)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**clumping**](#variable-clumping)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cosine**](#variable-cosine)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**coszen**](#variable-coszen)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fclai**](#variable-fclai)  <br> |
|  real(rk), dimension([**z**](canopy__var3din__mod_8F90.md#variable-z)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fractional**](#variable-fractional)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**index**](#variable-index)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lai**](#variable-lai)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**leaf**](#variable-leaf)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**model**](#variable-model)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**of**](#variable-of)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**shapes**](#variable-shapes)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**solar**](#variable-solar)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**the**](#variable-the)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**total**](#variable-total)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zenith**](#variable-zenith)  <br> |












































## Detailed Description


This subroutine computes the photolysis correction factors within a forest canopy using the exponential attenuation model from Makar et al. (2017). The calculation accounts for:
* Fractional cumulative leaf area index (FCLAI) profiles
* Total leaf area index (LAI)
* Clumping index to account for non-random leaf distribution
* Solar zenith angle effects on light penetration




The photolysis correction factor represents the fraction of photolysis rates relative to above-canopy conditions at each height within the canopy.




**Parameters:**


* `FCLAI` Fractional cumulative LAI shapes of plant surface distribution (nondimensional) 
* `LAI` Model input total Leaf Area Index (m²/m²) 
* `CLU` Model input Clumping Index (nondimensional) 
* `COSZEN` Model input Cosine Solar Zenith Angle (nondimensional) 
* `RJCF` Photolysis correction factor (nondimensional) 




    
## Public Attributes Documentation




### variable angle 

```Fortran
real(rk), intent(in) angle;
```




<hr>



### variable area 

```Fortran
real(rk), intent(in) area;
```




<hr>



### variable clu 

```Fortran
real(rk), intent(in) clu;
```




<hr>



### variable clumping 

```Fortran
real(rk), intent(in) clumping;
```




<hr>



### variable cosine 

```Fortran
real(rk), intent(in) cosine;
```




<hr>



### variable coszen 

```Fortran
real(rk), intent(in) coszen;
```




<hr>



### variable fclai 

```Fortran
real(rk), dimension(:), intent(in) fclai;
```




<hr>



### variable fractional 

```Fortran
real(rk), dimension (z), intent(in) fractional;
```




<hr>



### variable index 

```Fortran
real(rk), intent(in) index;
```




<hr>



### variable input 

```Fortran
real(rk), intent(in) input;
```




<hr>



### variable lai 

```Fortran
real(rk), intent(in) lai;
```




<hr>



### variable leaf 

```Fortran
real(rk), intent(in) leaf;
```




<hr>



### variable model 

```Fortran
real(rk), intent(in) model;
```




<hr>



### variable of 

```Fortran
real(rk), intent(in) of;
```




<hr>



### variable shapes 

```Fortran
real(rk), intent(in) shapes;
```




<hr>



### variable solar 

```Fortran
real(rk), intent(in) solar;
```




<hr>



### variable the 

```Fortran
real(rk), intent(in) the;
```




<hr>



### variable total 

```Fortran
real(rk), intent(in) total;
```




<hr>



### variable zenith 

```Fortran
real(rk), intent(in) zenith;
```




<hr>

------------------------------


