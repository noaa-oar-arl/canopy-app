

# Group write\_txt



[**Modules**](modules.md) **>** [**write\_txt**](group__write__txt.md)



_Routines for writing text format output files._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_write\_txt**](#function-canopy_write_txt) (OUTPREFX OUTPREFX, TIMENOW TIMENOW) <br>_Write canopy outputs to text files._  |




























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


