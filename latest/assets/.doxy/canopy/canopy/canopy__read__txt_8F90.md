

# File canopy\_read\_txt.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_read\_txt.F90**](canopy__read__txt_8F90.md)

[Go to the source code of this file](canopy__read__txt_8F90_source.md)

_Text Input File Reader Subroutine._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_read\_txt**](#function-canopy_read_txt) (INFILE INFILE, INFILE2 INFILE2) <br>_Read canopy met/sfc inputs from text files._  |




























## Detailed Description


This file contains the subroutine for reading meteorological and surface data from text format input files. It handles both main variable files and supplementary canopy profile files.




**Author:**

Patrick C. Campbell 




**Date:**

October 2022 





    
## Public Functions Documentation




### function canopy\_read\_txt 

_Read canopy met/sfc inputs from text files._ 
```Fortran
subroutine canopy_read_txt (
    INFILE INFILE,
    INFILE2 INFILE2
) 
```



This subroutine reads meteorological and surface input data from text files. It can handle both single variable text files and combined variable plus canopy profile text files depending on the var3d\_opt setting.




**Parameters:**


* `INFILE` Primary input text file path 
* `INFILE2` Secondary canopy profile input text file path 
* `infile` main IO text reader/writer 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_read_txt.F90`

