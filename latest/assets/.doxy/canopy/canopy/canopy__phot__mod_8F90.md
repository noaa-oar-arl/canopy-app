

# File canopy\_phot\_mod.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_phot\_mod.F90**](canopy__phot__mod_8F90.md)

[Go to the source code of this file](canopy__phot__mod_8F90_source.md)

_Canopy Photolysis Module._ [More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**canopy\_phot\_mod**](namespacecanopy__phot__mod.md) <br> |








## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**angle**](#variable-angle)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**area**](#variable-area)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**clu**](#variable-clu)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**clumping**](#variable-clumping)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**correction**](#variable-correction)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cosine**](#variable-cosine)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**coszen**](#variable-coszen)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**factor**](#variable-factor)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fclai**](#variable-fclai)  <br> |
|  real(rk), dimension([**z**](canopy__var3din__mod_8F90.md#variable-z)), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fractional**](#variable-fractional)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**index**](#variable-index)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lai**](#variable-lai)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**leaf**](#variable-leaf)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**model**](#variable-model)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**of**](#variable-of)  <br> |
|  real(rk), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**photolysis**](#variable-photolysis)  <br> |
|  real(rk), dimension(:), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**rjcf**](#variable-rjcf)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**shapes**](#variable-shapes)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**solar**](#variable-solar)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**the**](#variable-the)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**total**](#variable-total)  <br> |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zenith**](#variable-zenith)  <br> |












































## Detailed Description


This module contains subroutines for calculating photolysis attenuation within forest canopies based on the algorithms described in Makar et al. (2017). The module computes how photolysis rates are reduced within the canopy due to shading by leaves and branches.




**Author:**

Patrick C. Campbell 




**Date:**

June 2022


\references Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone. Nat Commun 8, 15243 (2017). [https://doi.org/10.1038/ncomms15243](https://doi.org/10.1038/ncomms15243) 


    
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



### variable correction 

```Fortran
real(rk), intent(out) correction;
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



### variable factor 

```Fortran
real(rk), intent(out) factor;
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



### variable photolysis 

```Fortran
real(rk), intent(out) photolysis;
```




<hr>



### variable rjcf 

```Fortran
real(rk), dimension(:), intent(out) rjcf;
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
The documentation for this class was generated from the following file `src/canopy_phot_mod.F90`

