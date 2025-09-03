

# Group canopy\_app



[**Modules**](modules.md) **>** [**canopy\_app**](group__canopy__app.md)



_Main program for the canopy model application._ 












## Modules

| Type | Name |
| ---: | :--- |
| module | [**Main Application Workflow**](group__main__workflow.md) <br>_Primary computational workflow of the canopy application._  |


























## Public Functions

| Type | Name |
| ---: | :--- |
|  program | [**canopy\_app**](#function-canopy_app) () <br>_Main Canopy Application Program._  |




























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


