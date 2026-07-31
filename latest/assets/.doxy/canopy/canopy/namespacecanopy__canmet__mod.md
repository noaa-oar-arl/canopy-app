

# Namespace canopy\_canmet\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_canmet\_mod**](namespacecanopy__canmet__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk) | [**canfracref**](#variable-canfracref)  <br>_Reference canopy fraction._  |
|  real(rk) | [**cluref**](#variable-cluref)  <br>_Reference clumping index._  |
|  real(rk) | [**cszref**](#variable-cszref)  <br>_Reference cosine of zenith angle._  |
|  real(rk) | [**dswrfref**](#variable-dswrfref)  <br>_Reference downward shortwave radiation._  |
|  real(rk) | [**frpref**](#variable-frpref)  <br>_Reference fire radiative power._  |
|  real(rk) | [**hcmref**](#variable-hcmref)  <br>_Reference canopy height._  |
|  real(rk) | [**hgtref**](#variable-hgtref)  <br>_Reference height above canopy._  |
|  real(rk) | [**hpblref**](#variable-hpblref)  <br>_Reference planetary boundary layer height._  |
|  real(rk) | [**icec\_averef**](#variable-icec_averef)  <br>_Reference average ground or water ice cover._  |
|  real(rk) | [**lairef**](#variable-lairef)  <br>_Reference leaf area index._  |
|  real(rk) | [**latref**](#variable-latref)  <br>_Met/Sfc variable reassignment names above reference conditions from the model._  |
|  real(rk), dimension(:), allocatable | [**lev\_arr**](#variable-lev_arr)  <br> |
|  real(rk), dimension(:), allocatable | [**levref**](#variable-levref)  <br>_Reference vertical levels array._  |
|  real(rk) | [**lonref**](#variable-lonref)  <br>_Reference longitude._  |
|  real(rk) | [**molref**](#variable-molref)  <br>_Reference Monin-Obukhov length._  |
|  real(rk) | [**ozone\_w126ref**](#variable-ozone_w126ref)  <br>_Reference ozone W126 values._  |
|  real(rk), dimension(:), allocatable | [**pavd\_arr**](#variable-pavd_arr)  <br> |
|  real(rk), dimension(:), allocatable | [**pavdref**](#variable-pavdref)  <br>_Reference plant area volume density array._  |
|  real(rk) | [**prate\_averef**](#variable-prate_averef)  <br>_Reference precipitation rate._  |
|  real(rk) | [**pressfcref**](#variable-pressfcref)  <br>_Reference surface pressure._  |
|  real(rk) | [**shtflref**](#variable-shtflref)  <br>_Reference surface sensible heat flux._  |
|  real(rk) | [**snowc\_averef**](#variable-snowc_averef)  <br>_Reference average ground snow cover._  |
|  real(rk) | [**soilt1ref**](#variable-soilt1ref)  <br>_Reference soil temperature level 1._  |
|  real(rk) | [**soilt2ref**](#variable-soilt2ref)  <br>_Reference soil temperature level 2._  |
|  real(rk) | [**soilt3ref**](#variable-soilt3ref)  <br>_Reference soil temperature level 3._  |
|  real(rk) | [**soilt4ref**](#variable-soilt4ref)  <br>_Reference soil temperature level 4._  |
|  real(rk) | [**soilw1ref**](#variable-soilw1ref)  <br>_Reference soil moisture layer 1._  |
|  real(rk) | [**soilw2ref**](#variable-soilw2ref)  <br>_Reference soil moisture layer 2._  |
|  real(rk) | [**soilw3ref**](#variable-soilw3ref)  <br>_Reference soil moisture layer 3._  |
|  real(rk) | [**soilw4ref**](#variable-soilw4ref)  <br>_Reference soil moisture layer 4._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**sotypref**](#variable-sotypref)  <br>_Reference soil type._  |
|  real(rk) | [**spfh2mref**](#variable-spfh2mref)  <br>_Reference 2-meter specific humidity._  |
|  real(rk) | [**tmp2mref**](#variable-tmp2mref)  <br>_Reference 2-meter temperature._  |
|  real(rk) | [**tmp\_hyblev1ref**](#variable-tmp_hyblev1ref)  <br>_Reference first model layer air temperature._  |
|  real(rk) | [**tmpsfcref**](#variable-tmpsfcref)  <br>_Reference surface temperature._  |
|  real(rk) | [**ubzref**](#variable-ubzref)  <br>_Reference bulk wind speed._  |
|  real(rk) | [**uref**](#variable-uref)  <br>_Reference U wind speed._  |
|  real(rk) | [**ustref**](#variable-ustref)  <br>_Reference friction velocity._  |
|  [**type**](canopy__bioparm__mod_8F90.md#variable-type)([**variable\_type**](namespacecanopy__canmet__mod.md#none-variable_type)), dimension(:), allocatable | [**variables**](#variable-variables)  <br>_Allocated array for 1D meteorological variables._  |
|  [**type**](canopy__bioparm__mod_8F90.md#variable-type)([**variable\_type\_1d**](namespacecanopy__canmet__mod.md#none-variable_type_1d)), dimension(:), allocatable | [**variables\_1d**](#variable-variables_1d)  <br>_Allocated array for 1D vertical variables._  |
|  [**type**](canopy__bioparm__mod_8F90.md#variable-type)([**variable\_type**](namespacecanopy__canmet__mod.md#none-variable_type)), dimension(:, :), allocatable | [**variables\_2d**](#variable-variables_2d)  <br> |
|  [**type**](canopy__bioparm__mod_8F90.md#variable-type)([**variable\_type\_3d**](namespacecanopy__canmet__mod.md#none-variable_type_3d)), dimension(:, :, :), allocatable | [**variables\_3d**](#variable-variables_3d)  <br>_Allocated array for 3D PAVD variables._  |
|  [**type**](canopy__bioparm__mod_8F90.md#variable-type)([**variable\_type\_can**](namespacecanopy__canmet__mod.md#none-variable_type_can)), dimension(:), allocatable | [**variables\_can**](#variable-variables_can)  <br>_Allocated array for canopy profile variables._  |
|  real(rk) | [**vbdsf\_averef**](#variable-vbdsf_averef)  <br>_Downward shortwave radiation._  |
|  real(rk) | [**vddsf\_averef**](#variable-vddsf_averef)  <br>_Downward shortwave radiation._  |
|  real(rk) | [**vref**](#variable-vref)  <br>_Reference V wind speed._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**vtyperef**](#variable-vtyperef)  <br>_Reference vegetation type._  |
|  real(rk) | [**wiltref**](#variable-wiltref)  <br>_Reference wilting point._  |
|  real(rk) | [**z0ref**](#variable-z0ref)  <br>_Reference surface roughness length._  |












































## Public Attributes Documentation




### variable canfracref 

_Reference canopy fraction._ 
```Fortran
real(rk) canopy_canmet_mod::canfracref;
```



Input canopy fraction of grid cell 

**Parameters:**


* `units` dimensionless fraction 




        

<hr>



### variable cluref 

_Reference clumping index._ 
```Fortran
real(rk) canopy_canmet_mod::cluref;
```



Input canopy clumping index 

**Parameters:**


* `units` dimensionless 




        

<hr>



### variable cszref 

_Reference cosine of zenith angle._ 
```Fortran
real(rk) canopy_canmet_mod::cszref;
```



Input cosine of zenith angle 

**Parameters:**


* `units` dimensionless 




        

<hr>



### variable dswrfref 

_Reference downward shortwave radiation._ 
```Fortran
real(rk) canopy_canmet_mod::dswrfref;
```



Instantaneous downward shortwave radiation 

**Parameters:**


* `units` W/m² 




        

<hr>



### variable frpref 

_Reference fire radiative power._ 
```Fortran
real(rk) canopy_canmet_mod::frpref;
```



Input fire radiative power 

**Parameters:**


* `units` MW 




        

<hr>



### variable hcmref 

_Reference canopy height._ 
```Fortran
real(rk) canopy_canmet_mod::hcmref;
```



Input canopy height 

**Parameters:**


* `units` meters (m) 




        

<hr>



### variable hgtref 

_Reference height above canopy._ 
```Fortran
real(rk) canopy_canmet_mod::hgtref;
```



Reference height above the canopy 

**Parameters:**


* `units` meters (m) 




        

<hr>



### variable hpblref 

_Reference planetary boundary layer height._ 
```Fortran
real(rk) canopy_canmet_mod::hpblref;
```



Height of planetary boundary layer 

**Parameters:**


* `units` meters (m) 




        

<hr>



### variable icec\_averef 

_Reference average ground or water ice cover._ 
```Fortran
real(rk) canopy_canmet_mod::icec_averef;
```



Average fraction ground or water ice cover 

**Parameters:**


* `units` dimensionless fraction 




        

<hr>



### variable lairef 

_Reference leaf area index._ 
```Fortran
real(rk) canopy_canmet_mod::lairef;
```



Input leaf area index 

**Parameters:**


* `units` m²/m² 




        

<hr>



### variable latref 

_Met/Sfc variable reassignment names above reference conditions from the model._ 
```Fortran
real(rk) canopy_canmet_mod::latref;
```



Reference meteorological and surface variables used in canopy calculations


Reference latitude


Latitude of cell/point 

**Parameters:**


* `units` degrees 




        

<hr>



### variable lev\_arr 

```Fortran
real(rk), dimension  ( : ), allocatable canopy_canmet_mod::lev_arr;
```




<hr>



### variable levref 

_Reference vertical levels array._ 
```Fortran
real(rk), dimension ( : ), allocatable canopy_canmet_mod::levref;
```



Reference vertical levels with 3D input data 

**Parameters:**


* `units` meters (m) 




        

<hr>



### variable lonref 

_Reference longitude._ 
```Fortran
real(rk) canopy_canmet_mod::lonref;
```



Longitude of cell/point 

**Parameters:**


* `units` degrees 




        

<hr>



### variable molref 

_Reference Monin-Obukhov length._ 
```Fortran
real(rk) canopy_canmet_mod::molref;
```



Input Monin-Obukhov Length 

**Parameters:**


* `units` meters (m) 




        

<hr>



### variable ozone\_w126ref 

_Reference ozone W126 values._ 
```Fortran
real(rk) canopy_canmet_mod::ozone_w126ref;
```



Ozone W126 values 

**Parameters:**


* `units` ppm-hours 




        

<hr>



### variable pavd\_arr 

```Fortran
real(rk), dimension ( : ), allocatable canopy_canmet_mod::pavd_arr;
```




<hr>



### variable pavdref 

_Reference plant area volume density array._ 
```Fortran
real(rk), dimension ( : ), allocatable canopy_canmet_mod::pavdref;
```



Plant area volume density 

**Parameters:**


* `units` m²/m³ 




        

<hr>



### variable prate\_averef 

_Reference precipitation rate._ 
```Fortran
real(rk) canopy_canmet_mod::prate_averef;
```



Mass precipitation rate 

**Parameters:**


* `units` kg/m²/s 




        

<hr>



### variable pressfcref 

_Reference surface pressure._ 
```Fortran
real(rk) canopy_canmet_mod::pressfcref;
```



Surface pressure 

**Parameters:**


* `units` hPa 




        

<hr>



### variable shtflref 

_Reference surface sensible heat flux._ 
```Fortran
real(rk) canopy_canmet_mod::shtflref;
```



Instantaneous surface sensible heat net flux 

**Parameters:**


* `units` W/m² 




        

<hr>



### variable snowc\_averef 

_Reference average ground snow cover._ 
```Fortran
real(rk) canopy_canmet_mod::snowc_averef;
```



Average percent ground snow cover 

**Parameters:**


* `units` percent (%) 




        

<hr>



### variable soilt1ref 

_Reference soil temperature level 1._ 
```Fortran
real(rk) canopy_canmet_mod::soilt1ref;
```



Soil temperature level 1 

**Parameters:**


* `units` K 




        

<hr>



### variable soilt2ref 

_Reference soil temperature level 2._ 
```Fortran
real(rk) canopy_canmet_mod::soilt2ref;
```



Soil temperature level 2 

**Parameters:**


* `units` K 




        

<hr>



### variable soilt3ref 

_Reference soil temperature level 3._ 
```Fortran
real(rk) canopy_canmet_mod::soilt3ref;
```



Soil temperature level 3 

**Parameters:**


* `units` K 




        

<hr>



### variable soilt4ref 

_Reference soil temperature level 4._ 
```Fortran
real(rk) canopy_canmet_mod::soilt4ref;
```



Soil temperature level 4 

**Parameters:**


* `units` K 




        

<hr>



### variable soilw1ref 

_Reference soil moisture layer 1._ 
```Fortran
real(rk) canopy_canmet_mod::soilw1ref;
```



Volumetric soil moisture layer 1 

**Parameters:**


* `units` m³/m³ 




        

<hr>



### variable soilw2ref 

_Reference soil moisture layer 2._ 
```Fortran
real(rk) canopy_canmet_mod::soilw2ref;
```



Volumetric soil moisture layer 2 

**Parameters:**


* `units` m³/m³ 




        

<hr>



### variable soilw3ref 

_Reference soil moisture layer 3._ 
```Fortran
real(rk) canopy_canmet_mod::soilw3ref;
```



Volumetric soil moisture layer 3 

**Parameters:**


* `units` m³/m³ 




        

<hr>



### variable soilw4ref 

_Reference soil moisture layer 4._ 
```Fortran
real(rk) canopy_canmet_mod::soilw4ref;
```



Volumetric soil moisture layer 4 

**Parameters:**


* `units` m³/m³ 




        

<hr>



### variable sotypref 

_Reference soil type._ 
```Fortran
integer canopy_canmet_mod::sotypref;
```



Soil type 


        

<hr>



### variable spfh2mref 

_Reference 2-meter specific humidity._ 
```Fortran
real(rk) canopy_canmet_mod::spfh2mref;
```



2-meter specific humidity 

**Parameters:**


* `units` kg/kg 




        

<hr>



### variable tmp2mref 

_Reference 2-meter temperature._ 
```Fortran
real(rk) canopy_canmet_mod::tmp2mref;
```



2-meter temperature 

**Parameters:**


* `units` K 




        

<hr>



### variable tmp\_hyblev1ref 

_Reference first model layer air temperature._ 
```Fortran
real(rk) canopy_canmet_mod::tmp_hyblev1ref;
```



1st model layer air temperature above ground 

**Parameters:**


* `units` K 




        

<hr>



### variable tmpsfcref 

_Reference surface temperature._ 
```Fortran
real(rk) canopy_canmet_mod::tmpsfcref;
```



Surface temperature 

**Parameters:**


* `units` K 




        

<hr>



### variable ubzref 

_Reference bulk wind speed._ 
```Fortran
real(rk) canopy_canmet_mod::ubzref;
```



Input above canopy/reference 10-m model wind speed 

**Parameters:**


* `units` m/s 




        

<hr>



### variable uref 

_Reference U wind speed._ 
```Fortran
real(rk) canopy_canmet_mod::uref;
```



Input above canopy/reference 10-m U wind speed 

**Parameters:**


* `units` m/s 




        

<hr>



### variable ustref 

_Reference friction velocity._ 
```Fortran
real(rk) canopy_canmet_mod::ustref;
```



Input friction velocity 

**Parameters:**


* `units` m/s 




        

<hr>



### variable variables 

_Allocated array for 1D meteorological variables._ 
```Fortran
type(variable_type), dimension( : ), allocatable canopy_canmet_mod::variables;
```



Array for storing 1D meteorological input variables 


        

<hr>



### variable variables\_1d 

_Allocated array for 1D vertical variables._ 
```Fortran
type(variable_type_1d), dimension( : ), allocatable canopy_canmet_mod::variables_1d;
```



Array for storing 1D vertical level variables 


        

<hr>



### variable variables\_2d 

```Fortran
type(variable_type), dimension( : , :), allocatable canopy_canmet_mod::variables_2d;
```




<hr>



### variable variables\_3d 

_Allocated array for 3D PAVD variables._ 
```Fortran
type(variable_type_3d), dimension( : , : , :), allocatable canopy_canmet_mod::variables_3d;
```



Array for storing 3D PAVD variables 


        

<hr>



### variable variables\_can 

_Allocated array for canopy profile variables._ 
```Fortran
type(variable_type_can), dimension( : ), allocatable canopy_canmet_mod::variables_can;
```



Array for storing canopy profile variables 


        

<hr>



### variable vbdsf\_averef 

_Downward shortwave radiation._ 
```Fortran
real(rk) canopy_canmet_mod::vbdsf_averef;
```



Average downward visible beam radiation 

**Parameters:**


* `units` W/m² 




        

<hr>



### variable vddsf\_averef 

_Downward shortwave radiation._ 
```Fortran
real(rk) canopy_canmet_mod::vddsf_averef;
```



Average downward visible diffuse radiation 

**Parameters:**


* `units` W/m² 




        

<hr>



### variable vref 

_Reference V wind speed._ 
```Fortran
real(rk) canopy_canmet_mod::vref;
```



Input above canopy/reference 10-m V wind speed 

**Parameters:**


* `units` m/s 




        

<hr>



### variable vtyperef 

_Reference vegetation type._ 
```Fortran
integer canopy_canmet_mod::vtyperef;
```



Input vegetation type (VIIRS) 


        

<hr>



### variable wiltref 

_Reference wilting point._ 
```Fortran
real(rk) canopy_canmet_mod::wiltref;
```



Wilting point 

**Parameters:**


* `units` proportion 




        

<hr>



### variable z0ref 

_Reference surface roughness length._ 
```Fortran
real(rk) canopy_canmet_mod::z0ref;
```



Input total/surface roughness length 

**Parameters:**


* `units` meters (m) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_canmet_mod.F90`

