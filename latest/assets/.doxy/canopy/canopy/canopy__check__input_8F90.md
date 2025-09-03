

# File canopy\_check\_input.F90



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_check\_input.F90**](canopy__check__input_8F90.md)

[Go to the source code of this file](canopy__check__input_8F90_source.md)

_Input File Validation and Reading Subroutine._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**canopy**](#variable-canopy)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**extension**](#variable-extension)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**file**](#variable-file)  <br> |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**infile**](#variable-infile)  <br> |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**infile2**](#variable-infile2)  <br> |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**input**](#variable-input)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**of**](#variable-of)  <br> |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**path**](#variable-path)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**position**](#variable-position)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**ppos**](#variable-ppos)  <br> |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**primary**](#variable-primary)  <br> |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**profile**](#variable-profile)  <br> |
|  character(len= \*), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**secondary**](#variable-secondary)  <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**separator**](#variable-separator)  <br> |












































## Detailed Description


This file contains the subroutine that checks and reads canopy input files. It supports both text and NetCDF input formats and validates the file format against the user-specified input format option.




**Author:**

Patrick C. Campbell 




**Date:**

December 2022 





    
## Public Attributes Documentation




### variable canopy 

```Fortran
character(len=*), intent(in) canopy;
```




<hr>



### variable extension 

```Fortran
integer extension;
```




<hr>



### variable file 

```Fortran
character(len=*), intent(in) file;
```




<hr>



### variable infile 

```Fortran
character(len=*), intent(in) infile;
```




<hr>



### variable infile2 

```Fortran
character(len=*), intent(in) infile2;
```




<hr>



### variable input 

```Fortran
character(len=*), intent(in) input;
```




<hr>



### variable of 

```Fortran
integer of;
```




<hr>



### variable path 

```Fortran
character(len=*), intent(in) path;
```




<hr>



### variable position 

```Fortran
integer position;
```




<hr>



### variable ppos 

```Fortran
integer ppos;
```




<hr>



### variable primary 

```Fortran
character(len=*), intent(in) primary;
```




<hr>



### variable profile 

```Fortran
character(len=*), intent(in) profile;
```




<hr>



### variable secondary 

```Fortran
character(len=*), intent(in) secondary;
```




<hr>



### variable separator 

```Fortran
integer separator;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_check_input.F90`

