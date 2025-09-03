

# Namespace canopy\_coord\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_coord\_mod**](namespacecanopy__coord__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**modlays**](#variable-modlays)  <br>_Number of total above and below canopy model layers._  |
|  real(rk) | [**modres**](#variable-modres)  <br>_Model above and below canopy vertical resolution [m]._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**nlat**](#variable-nlat)  <br>_Length of latitude coordinate (number of latitude points)_  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**nlon**](#variable-nlon)  <br>_Length of longitude coordinate (number of longitude points)_  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**ntime**](#variable-ntime)  <br>_Number of model timesteps._  |
|  character(len=24) | [**time\_end**](#variable-time_end)  <br>_Simulation end time in YYYY-MM-DD-HH:MM:SS.SSSS format._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**time\_intvl**](#variable-time_intvl)  <br>_Time interval for input/output [seconds]._  |
|  character(len=24) | [**time\_start**](#variable-time_start)  <br>_Simulation start time in YYYY-MM-DD-HH:MM:SS.SSSS format._  |












































## Public Attributes Documentation




### variable modlays 

_Number of total above and below canopy model layers._ 
```Fortran
integer canopy_coord_mod::modlays;
```




<hr>



### variable modres 

_Model above and below canopy vertical resolution [m]._ 
```Fortran
real(rk) canopy_coord_mod::modres;
```




<hr>



### variable nlat 

_Length of latitude coordinate (number of latitude points)_ 
```Fortran
integer canopy_coord_mod::nlat;
```




<hr>



### variable nlon 

_Length of longitude coordinate (number of longitude points)_ 
```Fortran
integer canopy_coord_mod::nlon;
```




<hr>



### variable ntime 

_Number of model timesteps._ 
```Fortran
integer canopy_coord_mod::ntime;
```




<hr>



### variable time\_end 

_Simulation end time in YYYY-MM-DD-HH:MM:SS.SSSS format._ 
```Fortran
character(len=24) canopy_coord_mod::time_end;
```




<hr>



### variable time\_intvl 

_Time interval for input/output [seconds]._ 
```Fortran
integer canopy_coord_mod::time_intvl;
```




<hr>



### variable time\_start 

_Simulation start time in YYYY-MM-DD-HH:MM:SS.SSSS format._ 
```Fortran
character(len=24) canopy_coord_mod::time_start;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_coord_mod.F90`

