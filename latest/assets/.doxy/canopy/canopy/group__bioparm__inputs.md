

# Group bioparm\_inputs



[**Modules**](modules.md) **>** [**bioparm\_inputs**](group__bioparm__inputs.md)



_Get biogenic emission factors and parameters from MEGAN2.1._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), dimension(default=0/viirs), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**al**](#variable-al)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**biogenic**](#variable-biogenic)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**cell**](#variable-cell)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**dominant**](#variable-dominant)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**emi\_ind**](#variable-emi_ind)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**emissions**](#variable-emissions)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**et**](#variable-et)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**for**](#variable-for)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**from**](#variable-from)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**grid**](#variable-grid)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**index**](#variable-index)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  integer, intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**integer**](#variable-integer)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lu**](#variable-lu)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**lu\_opt**](#variable-lu_opt)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**mapped**](#variable-mapped)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**massman**](#variable-massman)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**model**](#variable-model)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**to**](#variable-to)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**type**](#variable-type)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vegetation**](#variable-vegetation)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**vtype**](#variable-vtype)  <br> |












































## Detailed Description


This subroutine retrieves biogenic emission factors and parameters based on the MEGAN2.1 framework. It provides plant-dependent emission capacities for various biogenic volatile organic compounds including:
* Isoprene
* Myrcene
* Sabinene
* Limonene
* 3-Carene
* T-β-Ocimene
* β-Pinene
* α-Pinene
* 2-Methyl-3-buten-2-ol (MBO)
* Methanol
* Acetone
* Other monoterpenes and sesquiterpenes




The parameters are provided for different vegetation functional types and include emission factors, light-dependent fractions, temperature coefficients, leaf age factors, and stress response parameters.




**Parameters:**


* `EMI_IND` Input biogenic emissions index 
* `LU_OPT` Land use type option (0=VIIRS, 1=MODIS) 
* `VTYPE` Grid cell dominant vegetation type 
* `EF` Mapped emission factor (μg/m²/hr) 
* `LDF` Light-dependent fraction 
* `BETA` Empirical coefficient for temperature dependence of light-independent fraction 
* `CT1` Activation energy (kJ/mol) 
* `CEO` Empirical coefficient 
* `ANEW` Empirical factor for new foliage 
* `AGRO` Empirical factor for growing foliage 
* `AMAT` Empirical factor for mature foliage 
* `AOLD` Empirical factor for old/senescing foliage 
* `ROOTA` Coefficient A for PFT dependent cumulative root depth fraction (m⁻¹) 
* `ROOTB` Coefficient B for PFT dependent cumulative root depth fraction (m⁻¹) 
* `CAQ` Coefficient for poor air quality stress 
* `TAQ` Threshold for poor air quality stress (ppm-hours) 
* `DTAQ` Delta threshold for poor air quality stress (ppm-hours) 
* `CHT` Coefficient for high temperature stress 
* `THT` Threshold for high temperature stress (K) 
* `DTHT` Delta threshold for high temperature stress (K) 
* `CLT` Coefficient for low temperature stress 
* `TLT` Threshold for low temperature stress (K) 
* `DTLT` Delta threshold for low temperature stress (K) 
* `CHW` Coefficient for high wind stress 
* `THW` Threshold for high wind stress (m/s) 
* `DTHW` Delta threshold for high wind stress (m/s) 




    
## Public Attributes Documentation




### variable al 

```Fortran
integer, dimension (default = 0/viirs), intent(in) al;
```




<hr>



### variable biogenic 

```Fortran
integer, intent(in) biogenic;
```




<hr>



### variable cell 

```Fortran
integer, intent(in) cell;
```




<hr>



### variable dominant 

```Fortran
integer, intent(in) dominant;
```




<hr>



### variable emi\_ind 

```Fortran
integer, intent(in) emi_ind;
```




<hr>



### variable emissions 

```Fortran
integer, intent(in) emissions;
```




<hr>



### variable et 

```Fortran
integer, intent(in) et;
```




<hr>



### variable for 

```Fortran
integer, intent(in) for;
```




<hr>



### variable from 

```Fortran
integer, intent(in) from;
```




<hr>



### variable grid 

```Fortran
integer, intent(in) grid;
```




<hr>



### variable index 

```Fortran
integer, intent(in) index;
```




<hr>



### variable input 

```Fortran
integer, intent(in) input;
```




<hr>



### variable integer 

```Fortran
integer, intent(in) integer;
```




<hr>



### variable lu 

```Fortran
integer, intent(in) lu;
```




<hr>



### variable lu\_opt 

```Fortran
integer, intent(in) lu_opt;
```




<hr>



### variable mapped 

```Fortran
real(rk), intent(out) mapped;
```




<hr>



### variable massman 

```Fortran
integer, intent(in) massman;
```




<hr>



### variable model 

```Fortran
integer, intent(in) model;
```




<hr>



### variable to 

```Fortran
integer, intent(in) to;
```




<hr>



### variable type 

```Fortran
integer, intent(in) type;
```




<hr>



### variable vegetation 

```Fortran
integer, intent(in) vegetation;
```




<hr>



### variable vtype 

```Fortran
integer, intent(in) vtype;
```




<hr>

------------------------------


