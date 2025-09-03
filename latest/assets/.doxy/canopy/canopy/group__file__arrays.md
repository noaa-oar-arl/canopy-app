

# Group file\_arrays



[**Modules**](modules.md) **>** [**file\_arrays**](group__file__arrays.md)



_Arrays storing file paths for different input/output types._ 






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  character(len=256), dimension([**max\_mm**](namespacecanopy__files__mod.md#variable-max_mm)) | [**canopy\_files\_mod::file\_canvars**](#variable-file_canvars)  <br>_Array of canopy variable input file paths._  |
|  character(len=256), dimension(1) | [**canopy\_files\_mod::file\_out**](#variable-file_out)  <br>_Array of output file paths._  |
|  character(len=256), dimension([**max\_mm**](namespacecanopy__files__mod.md#variable-max_mm)) | [**canopy\_files\_mod::file\_vars**](#variable-file_vars)  <br>_Array of main variable input file paths._  |












































## Public Attributes Documentation




### variable file\_canvars 

_Array of canopy variable input file paths._ 
```Fortran
character(len=256), dimension ( max_mm ) canopy_files_mod::file_canvars;
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

------------------------------


