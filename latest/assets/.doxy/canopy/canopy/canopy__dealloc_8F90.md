

# File canopy\_dealloc.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_dealloc.F90**](canopy__dealloc_8F90.md)

[Go to the source code of this file](canopy__dealloc_8F90_source.md)

_Deallocation subroutine for canopy model arrays._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_dealloc**](#function-canopy_dealloc) () <br>_Deallocate arrays for canopy model._  |




























## Detailed Description


This subroutine deallocates all allocated arrays used in the canopy model to free memory at the end of model execution. 

**Author:**

P.C. Campbell 




**Date:**

October 2022 




**Version:**

1.0 





    
## Public Functions Documentation




### function canopy\_dealloc 

_Deallocate arrays for canopy model._ 
```Fortran
subroutine canopy_dealloc () 
```



Deallocates all allocated arrays including:
* Input variable arrays
* Canopy distribution arrays
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
The documentation for this class was generated from the following file `src/canopy_dealloc.F90`

