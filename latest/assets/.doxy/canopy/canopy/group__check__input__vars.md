

# Group check\_input\_vars



[**Modules**](modules.md) **>** [**check\_input\_vars**](group__check__input__vars.md)



_Check and read canopy input files (TXT or NETCDF)_ [More...](#detailed-description)






















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


This subroutine determines the input file format based on the file extension and validates it against the user-specified format option (infmt\_opt). It then calls the appropriate reading routine:
* For .txt files: calls [**canopy\_read\_txt()**](canopy__read__txt_8F90.md#function-canopy_read_txt)
* For .nc/.ncf/.nc4 files: calls [**canopy\_read\_ncf()**](namespacecanopy__ncf__io__mod.md#function-canopy_read_ncf)




The subroutine performs error checking to ensure the file format matches the namelist specification and exits with an error code if there are mismatches.




**Parameters:**


* `INFILE` Primary input file path 
* `INFILE2` Secondary canopy profile input file path 




    
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


