

# Namespace canopy\_files\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_files\_mod**](namespacecanopy__files__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**cdfid\_m**](#variable-cdfid_m)  <br>_NetCDF file identifier for main input files._  |
|  character(len=256), dimension([**max\_mm**](namespacecanopy__files__mod.md#variable-max_mm)) | [**file\_canvars**](#variable-file_canvars)  <br>_Array of canopy variable input file paths._  |
|  character(len= \*), parameter | [**file\_nml**](#variable-file_nml)   = `'[**input**](canopy__var3din__mod_8F90.md#variable-input)/[**namelist.canopy**](canopy__var3din__mod_8F90.md#variable-canopy)'`<br>_Path to the namelist configuration file._  |
|  character(len=256), dimension(1) | [**file\_out**](#variable-file_out)  <br>_Array of output file paths._  |
|  character(len=256), dimension([**max\_mm**](namespacecanopy__files__mod.md#variable-max_mm)) | [**file\_vars**](#variable-file_vars)  <br>_Array of main variable input file paths._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), parameter | [**iutnml**](#variable-iutnml)   = `8`<br>_FORTRAN unit number for namelist file._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), parameter | [**max\_mm**](#variable-max_mm)   = `10000`<br>_Maximum number of input files allowed._  |












































## Public Attributes Documentation




### variable cdfid\_m 

_NetCDF file identifier for main input files._ 
```Fortran
integer canopy_files_mod::cdfid_m;
```




<hr>



### variable file\_canvars 

_Array of canopy variable input file paths._ 
```Fortran
character(len=256), dimension ( max_mm ) canopy_files_mod::file_canvars;
```




<hr>



### variable file\_nml 

_Path to the namelist configuration file._ 
```Fortran
character(len=*), parameter canopy_files_mod::file_nml;
```




<hr>



### variable file\_out 

_Array of output file paths._ 
```Fortran
character(len=256), dimension     ( 1 ) canopy_files_mod::file_out;
```




<hr>



### variable file\_vars 

_Array of main variable input file paths._ 
```Fortran
character(len=256), dimension    ( max_mm ) canopy_files_mod::file_vars;
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
The documentation for this class was generated from the following file `src/canopy_files_mod.F90`

