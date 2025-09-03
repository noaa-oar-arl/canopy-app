

# Group read\_txt



[**Modules**](modules.md) **>** [**read\_txt**](group__read__txt.md)



_Routines for reading text format input files._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_read\_txt**](#function-canopy_read_txt) (INFILE INFILE, INFILE2 INFILE2) <br>_Read canopy met/sfc inputs from text files._  |




























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


