

# Namespace canopy\_ncf\_io\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_ncf\_io\_mod**](namespacecanopy__ncf__io__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**canopy\_close\_files**](#function-canopy_close_files) (OUTPREFX OUTPREFX) <br> |
|  subroutine, public | [**canopy\_outncf\_alloc**](#function-canopy_outncf_alloc) () <br> |
|  subroutine, public | [**canopy\_outncf\_init**](#function-canopy_outncf_init) () <br> |
|  subroutine, public | [**canopy\_outncfglobal**](#function-canopy_outncfglobal) (cdfid\_in cdfid\_in, fl fl) <br> |
|  subroutine, public | [**canopy\_read\_ncf**](#function-canopy_read_ncf) ([**infile**](canopy__check__input_8F90.md#variable-infile) infile) <br> |
|  subroutine, public | [**canopy\_write\_ncf**](#function-canopy_write_ncf) (OUTPREFX OUTPREFX) <br> |
|  subroutine | [**get\_var\_1d\_int\_cdf**](#function-get_var_1d_int_cdf) (cdfid cdfid, var var, idum1d idum1d, it it, rcode rcode) <br> |
|  subroutine | [**get\_var\_1d\_real\_cdf**](#function-get_var_1d_real_cdf) (cdfid cdfid, var var, dum1d dum1d, it it, rcode rcode) <br> |
|  subroutine | [**get\_var\_2d\_int\_cdf**](#function-get_var_2d_int_cdf) (cdfid cdfid, var var, idum2d idum2d, it it, rcode rcode) <br> |
|  subroutine | [**get\_var\_2d\_real\_cdf**](#function-get_var_2d_real_cdf) (cdfid cdfid, var var, dum2d dum2d, it it, rcode rcode) <br> |




























## Public Functions Documentation




### function canopy\_close\_files 

```Fortran
subroutine, public canopy_ncf_io_mod::canopy_close_files (
    OUTPREFX OUTPREFX
) 
```




<hr>



### function canopy\_outncf\_alloc 

```Fortran
subroutine, public canopy_ncf_io_mod::canopy_outncf_alloc () 
```




<hr>



### function canopy\_outncf\_init 

```Fortran
subroutine, public canopy_ncf_io_mod::canopy_outncf_init () 
```




<hr>



### function canopy\_outncfglobal 

```Fortran
subroutine, public canopy_ncf_io_mod::canopy_outncfglobal (
    cdfid_in cdfid_in,
    fl fl
) 
```




<hr>



### function canopy\_read\_ncf 

```Fortran
subroutine, public canopy_ncf_io_mod::canopy_read_ncf (
    infile infile
) 
```




<hr>



### function canopy\_write\_ncf 

```Fortran
subroutine, public canopy_ncf_io_mod::canopy_write_ncf (
    OUTPREFX OUTPREFX
) 
```




<hr>



### function get\_var\_1d\_int\_cdf 

```Fortran
subroutine canopy_ncf_io_mod::get_var_1d_int_cdf (
    cdfid cdfid,
    var var,
    idum1d idum1d,
    it it,
    rcode rcode
) 
```




<hr>



### function get\_var\_1d\_real\_cdf 

```Fortran
subroutine canopy_ncf_io_mod::get_var_1d_real_cdf (
    cdfid cdfid,
    var var,
    dum1d dum1d,
    it it,
    rcode rcode
) 
```




<hr>



### function get\_var\_2d\_int\_cdf 

```Fortran
subroutine canopy_ncf_io_mod::get_var_2d_int_cdf (
    cdfid cdfid,
    var var,
    idum2d idum2d,
    it it,
    rcode rcode
) 
```




<hr>



### function get\_var\_2d\_real\_cdf 

```Fortran
subroutine canopy_ncf_io_mod::get_var_2d_real_cdf (
    cdfid cdfid,
    var var,
    dum2d dum2d,
    it it,
    rcode rcode
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_ncf_io_mod.F90`

