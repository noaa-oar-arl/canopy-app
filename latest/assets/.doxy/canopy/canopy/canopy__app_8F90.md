

# File canopy\_app.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_app.F90**](canopy__app_8F90.md)

[Go to the source code of this file](canopy__app_8F90_source.md)

_Main Canopy Application Program._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  program | [**canopy\_app**](#function-canopy_app) () <br>_Main Canopy Application Program._  |




























## Detailed Description


This is the main application program that coordinates the canopy model calculations by orchestrating input reading, model initialization, time stepping, canopy calculations, and output writing.




**Author:**

Patrick C. Campbell 




**Date:**

June 2022




**Version:**


* Prototype: Patrick C. Campbell, 06/2022
* Revised: PCC (10/2022)
* 21 Aug 2023: Adding multiple timesteps (P.C. Campbell)
* Revised: 30 Nov 2023: Added supplementary canopy profile, file\_canvars (P.C. Campbell) 







    
## Public Functions Documentation




### function canopy\_app 

_Main Canopy Application Program._ 
```Fortran
program canopy_app () 
```



This program coordinates the entire canopy model workflow including:
* Reading user options from namelist
* Allocating and initializing variables
* Setting up input/output data structures
* Time stepping through the simulation period
* Calling main canopy calculations
* Writing output files
* Cleanup and deallocation 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_app.F90`

