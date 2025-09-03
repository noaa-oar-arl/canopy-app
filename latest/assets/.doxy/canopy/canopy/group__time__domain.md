

# Group time\_domain



[**Modules**](modules.md) **>** [**time\_domain**](group__time__domain.md)



_Variables defining the temporal domain of the simulation._ 






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**canopy\_coord\_mod::ntime**](#variable-ntime)  <br>_Number of model timesteps._  |
|  character(len=24) | [**canopy\_coord\_mod::time\_end**](#variable-time_end)  <br>_Simulation end time in YYYY-MM-DD-HH:MM:SS.SSSS format._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**canopy\_coord\_mod::time\_intvl**](#variable-time_intvl)  <br>_Time interval for input/output [seconds]._  |
|  character(len=24) | [**canopy\_coord\_mod::time\_start**](#variable-time_start)  <br>_Simulation start time in YYYY-MM-DD-HH:MM:SS.SSSS format._  |












































## Public Attributes Documentation




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


