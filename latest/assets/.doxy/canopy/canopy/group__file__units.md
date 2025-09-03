

# Group file\_units



[**Modules**](modules.md) **>** [**file\_units**](group__file__units.md)



_File unit numbers and array size limits._ 






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**canopy\_files\_mod::cdfid\_m**](#variable-cdfid_m)  <br>_NetCDF file identifier for main input files._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), parameter | [**canopy\_files\_mod::iutnml**](#variable-iutnml)   = `8`<br>_FORTRAN unit number for namelist file._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), parameter | [**canopy\_files\_mod::max\_mm**](#variable-max_mm)   = `10000`<br>_Maximum number of input files allowed._  |












































## Public Attributes Documentation




### variable cdfid\_m 

_NetCDF file identifier for main input files._ 
```Fortran
integer canopy_files_mod::cdfid_m;
```




<hr>



### variable iutnml 

_FORTRAN unit number for namelist file._ 
```Fortran
integer, parameter canopy_files_mod::iutnml;
```




<hr>



### variable max\_mm 

_Maximum number of input files allowed._ 
```Fortran
integer, parameter canopy_files_mod::max_mm;
```




<hr>

------------------------------


