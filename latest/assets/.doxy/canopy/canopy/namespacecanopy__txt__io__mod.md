

# Namespace canopy\_txt\_io\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_txt\_io\_mod**](namespacecanopy__txt__io__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**write\_txt**](#function-write_txt) (TXTPREFX TXTPREFX, TIMENOW TIMENOW) <br>_Write canopy model output to text file._  |




























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
The documentation for this class was generated from the following file `src/canopy_txt_io_mod.F90`

