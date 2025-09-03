

# Group txtio\_group



[**Modules**](modules.md) **>** [**txtio\_group**](group__txtio__group.md)



_Routines for reading and writing text format data files._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**i0**](#variable-i0)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**loc**](#variable-loc)  <br>_I/O status and location counter._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_txt\_io\_mod::write\_txt**](#function-write_txt) (TXTPREFX TXTPREFX, TIMENOW TIMENOW) <br>_Write canopy model output to text file._  |




























## Detailed Description


This group contains subroutines for reading meteorological and canopy variable text files and writing model output to text format. 


    
## Public Attributes Documentation




### variable i0 

```Fortran
integer i0;
```




<hr>



### variable loc 

_I/O status and location counter._ 
```Fortran
integer loc;
```




<hr>
## Public Functions Documentation




### function write\_txt 

_Write canopy model output to text file._ 
```Fortran
subroutine canopy_txt_io_mod::write_txt (
    TXTPREFX TXTPREFX,
    TIMENOW TIMENOW
) 
```



Writes comprehensive canopy model output variables to a text file with proper formatting and headers. Includes meteorological, canopy, and calculated variables for the specified time. 

**Parameters:**


* `TXTPREFX` Prefix for output text filename 
* `TIMENOW` Current simulation time string 



**Author:**

P.C. Campbell 




**Date:**

October 2022 





        

<hr>

------------------------------


