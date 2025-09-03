

# File canopy\_init.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_init.F90**](canopy__init_8F90.md)

[Go to the source code of this file](canopy__init_8F90_source.md)

_Initialization subroutine for canopy model arrays._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_init**](#function-canopy_init) () <br>_Initialize arrays for canopy model inputs and outputs._  |




























## Detailed Description


This subroutine initializes all arrays for canopy model inputs and outputs, setting them to appropriate fill values or zeros. 

**Author:**

P.C. Campbell 




**Date:**

October 2022 




**Version:**

1.0 





    
## Public Functions Documentation




### function canopy\_init 

_Initialize arrays for canopy model inputs and outputs._ 
```Fortran
subroutine canopy_init () 
```



Initializes all arrays for canopy model including:
* Canopy distribution variables
* Meteorological arrays
* Leaf temperature and PPFD arrays
* 24-hour and 240-hour temporary arrays
* Biogenic emission arrays
* Chemical species arrays 

**Author:**

P.C. Campbell 




**Date:**

October 2022 







        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_init.F90`

