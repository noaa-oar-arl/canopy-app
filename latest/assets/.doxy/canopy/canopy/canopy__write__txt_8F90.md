

# File canopy\_write\_txt.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_write\_txt.F90**](canopy__write__txt_8F90.md)

[Go to the source code of this file](canopy__write__txt_8F90_source.md)

_Text Output File Writer Subroutine._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_write\_txt**](#function-canopy_write_txt) (OUTPREFX OUTPREFX, TIMENOW TIMENOW) <br>_Write canopy outputs to text files._  |




























## Detailed Description


This file contains the subroutine for writing canopy model output to text format files. It serves as a wrapper that calls the main text writing routine with appropriate parameters.




**Author:**

Patrick C. Campbell 




**Date:**

October 2022 





    
## Public Functions Documentation




### function canopy\_write\_txt 

_Write canopy outputs to text files._ 
```Fortran
subroutine canopy_write_txt (
    OUTPREFX OUTPREFX,
    TIMENOW TIMENOW
) 
```



This subroutine writes canopy model calculation results to text format output files. It serves as a wrapper for the main text writing routine, passing the output file prefix and current time information.




**Parameters:**


* `OUTPREFX` Output file prefix string 
* `TIMENOW` Current time stamp string 
* `outprefx` main IO text reader/writer 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_write_txt.F90`

