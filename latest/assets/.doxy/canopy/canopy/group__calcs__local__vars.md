

# Group calcs\_local\_vars



[**Modules**](modules.md) **>** [**calcs\_local\_vars**](group__calcs__local__vars.md)



_Main canopy calculations dependent on canopy conditions._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk) | [**aerodynamic**](#variable-aerodynamic)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**and**](#variable-and)  <br> |
|  real(rk) | [**arrays**](#variable-arrays)  <br> |
|  real(rk) | [**averaging**](#variable-averaging)  <br> |
|  real(rk) | [**between**](#variable-between)  <br> |
|  real(rk) | [**biogenics**](#variable-biogenics)  <br> |
|  real(rk), dimension(saved [**from**](canopy__var3din__mod_8F90.md#variable-from) one [**timestep**](canopy__calcs_8F90.md#variable-timestep) [**to**](canopy__bioparm__mod_8F90.md#variable-to) [**the**](canopy__phot__mod_8F90.md#variable-the) next), save | [**cm2**](#variable-cm2)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**counters**](#variable-counters)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**current**](#variable-current)  <br> |
|  real(rk), save | [**currentlai**](#variable-currentlai)  <br> |
|  real(rk) | [**d**](#variable-d)  <br> |
|  real(rk) | [**days**](#variable-days)  <br> |
|  real(rk) | [**dep**](#variable-dep)  <br> |
|  real(rk) | [**dnewfrac**](#variable-dnewfrac)  <br> |
|  real(rk) | [**doldfrac**](#variable-doldfrac)  <br> |
|  real(rk) | [**dry**](#variable-dry)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**elapsed**](#variable-elapsed)  <br> |
|  real(rk) | [**for**](#variable-for)  <br> |
|  real(rk) | [**gas**](#variable-gas)  <br> |
|  real(rk) | [**historical**](#variable-historical)  <br> |
|  real(rk) | [**hnewfrac**](#variable-hnewfrac)  <br> |
|  real(rk) | [**holdfrac**](#variable-holdfrac)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**i**](#variable-i)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**in**](#variable-in)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**index**](#variable-index)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**int**](#variable-int)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**int\_nlaic**](#variable-int_nlaic)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**int\_nlaip**](#variable-int_nlaip)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**j**](#variable-j)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**k**](#variable-k)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**lai**](#variable-lai)  <br> |
|  real(rk), dimension(nlon \*nlat) | [**lat1d**](#variable-lat1d)  <br> |
|  real(rk), dimension(nlon, nlat) | [**lat2d**](#variable-lat2d)  <br> |
|  real(rk) | [**latitude**](#variable-latitude)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**loc**](#variable-loc)  <br>_I/O status and location counter._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**location**](#variable-location)  <br> |
|  real(rk), dimension(nlon \*nlat) | [**lon1d**](#variable-lon1d)  <br> |
|  real(rk), dimension(nlon, nlat) | [**lon2d**](#variable-lon2d)  <br> |
|  real(rk) | [**longitude**](#variable-longitude)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**loop**](#variable-loop)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**model**](#variable-model)  <br> |
|  real(rk) | [**nlaic**](#variable-nlaic)  <br> |
|  real(rk) | [**nlaip**](#variable-nlaip)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**number**](#variable-number)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**of**](#variable-of)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**past**](#variable-past)  <br> |
|  real(rk), save | [**pastlai**](#variable-pastlai)  <br> |
|  real(rk) | [**ra**](#variable-ra)  <br> |
|  real(rk) | [**resistance**](#variable-resistance)  <br> |
|  real(rk) | [**rib**](#variable-rib)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**the**](#variable-the)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**timestep**](#variable-timestep)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), save | [**timesteps**](#variable-timesteps)  <br> |
|  real(rk) | [**tsteplai**](#variable-tsteplai)  <br> |
|  real(rk) | [**variables**](#variable-variables)  <br> |












































## Detailed Description


This subroutine contains the main computational workflow for the canopy model. It processes meteorological inputs, performs canopy parameter calculations, and computes various canopy processes including:
* Grid distance calculations for 2D domains
* Reference height assignments
* Canopy geometry and morphology calculations
* Wind profile calculations within and above canopy
* Radiation transfer and photosynthesis computations
* Leaf temperature calculations for sun and shade conditions
* Biogenic emission calculations
* Dry deposition velocity computations
* Fire-related wind adjustment factors
* Photolysis attenuation factors
* Eddy diffusivity profiles




The subroutine handles both 1D point calculations and 2D gridded computations depending on the input format option specified in the namelist.




**Parameters:**


* `nn` Input time step index 
* `nn` aerosol dry deposition calculations 




    
## Public Attributes Documentation




### variable aerodynamic 

```Fortran
real(rk) aerodynamic;
```




<hr>



### variable and 

```Fortran
real(rk) and;
```




<hr>



### variable arrays 

```Fortran
real(rk) arrays;
```




<hr>



### variable averaging 

```Fortran
real(rk) averaging;
```




<hr>



### variable between 

```Fortran
real(rk) between;
```




<hr>



### variable biogenics 

```Fortran
real(rk) biogenics;
```




<hr>



### variable cm2 

```Fortran
real(rk), dimension  (saved from one timestep to the next), save cm2;
```




<hr>



### variable counters 

```Fortran
integer counters;
```




<hr>



### variable current 

```Fortran
real(rk), save current;
```




<hr>



### variable currentlai 

```Fortran
real(rk), save currentlai;
```




<hr>



### variable d 

```Fortran
real(rk) d;
```




<hr>



### variable days 

```Fortran
real(rk) days;
```




<hr>



### variable dep 

```Fortran
real(rk) dep;
```




<hr>



### variable dnewfrac 

```Fortran
real(rk) dnewfrac;
```




<hr>



### variable doldfrac 

```Fortran
real(rk) doldfrac;
```




<hr>



### variable dry 

```Fortran
real(rk) dry;
```




<hr>



### variable elapsed 

```Fortran
integer save elapsed;
```




<hr>



### variable for 

```Fortran
real(rk) for;
```




<hr>



### variable gas 

```Fortran
real(rk) gas;
```




<hr>



### variable historical 

```Fortran
real(rk) historical;
```




<hr>



### variable hnewfrac 

```Fortran
real(rk) hnewfrac;
```




<hr>



### variable holdfrac 

```Fortran
real(rk) holdfrac;
```




<hr>



### variable i 

```Fortran
integer i;
```




<hr>



### variable in 

```Fortran
integer save in;
```




<hr>



### variable index 

```Fortran
integer index;
```




<hr>



### variable int 

```Fortran
integer save int;
```




<hr>



### variable int\_nlaic 

```Fortran
integer int_nlaic;
```




<hr>



### variable int\_nlaip 

```Fortran
integer, save int_nlaip;
```




<hr>



### variable j 

```Fortran
integer j;
```




<hr>



### variable k 

```Fortran
integer k;
```




<hr>



### variable lai 

```Fortran
real(rk), save lai;
```




<hr>



### variable lat1d 

```Fortran
real(rk), dimension(nlon*nlat) lat1d;
```




<hr>



### variable lat2d 

```Fortran
real(rk), dimension(nlon,nlat) lat2d;
```




<hr>



### variable latitude 

```Fortran
real(rk) latitude;
```




<hr>



### variable loc 

_I/O status and location counter._ 
```Fortran
integer loc;
```




<hr>



### variable location 

```Fortran
integer location;
```




<hr>



### variable lon1d 

```Fortran
real(rk), dimension(nlon*nlat) lon1d;
```




<hr>



### variable lon2d 

```Fortran
real(rk), dimension(nlon,nlat) lon2d;
```




<hr>



### variable longitude 

```Fortran
real(rk) longitude;
```




<hr>



### variable loop 

```Fortran
integer loop;
```




<hr>



### variable model 

```Fortran
real(rk), save model;
```




<hr>



### variable nlaic 

```Fortran
real(rk) nlaic;
```




<hr>



### variable nlaip 

```Fortran
real(rk) nlaip;
```




<hr>



### variable number 

```Fortran
real(rk), save number;
```




<hr>



### variable of 

```Fortran
real(rk), save of;
```




<hr>



### variable past 

```Fortran
real(rk), save past;
```




<hr>



### variable pastlai 

```Fortran
real(rk), save pastlai;
```




<hr>



### variable ra 

```Fortran
real(rk) ra;
```




<hr>



### variable resistance 

```Fortran
real(rk) resistance;
```




<hr>



### variable rib 

```Fortran
real(rk) rib;
```




<hr>



### variable the 

```Fortran
real(rk), save the;
```




<hr>



### variable timestep 

```Fortran
integer save timestep;
```




<hr>



### variable timesteps 

```Fortran
integer save timesteps;
```




<hr>



### variable tsteplai 

```Fortran
real(rk) tsteplai;
```




<hr>



### variable variables 

```Fortran
real(rk) variables;
```




<hr>

------------------------------


