

# Group date\_group



[**Modules**](modules.md) **>** [**date\_group**](group__date__group.md)



_Routines for date and time manipulation and conversion._ [More...](#detailed-description)






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_date\_mod::geth\_idts**](#function-geth_idts) (newdate newdate, olddate olddate, idts idts) <br>_Calculate time difference between two dates._  |
|  subroutine | [**canopy\_date\_mod::geth\_newdate**](#function-geth_newdate) (ndate ndate, odate odate, idt idt) <br> |
|  subroutine | [**canopy\_date\_mod::getsdt**](#function-getsdt) (hdate hdate, sdate sdate, stime stime) <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) [**function**](canopy__var3din__mod_8F90.md#variable-function) | [**canopy\_date\_mod::julian**](#function-julian) (year year, mnth mnth, mday mday) <br> |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) [**function**](canopy__var3din__mod_8F90.md#variable-function) | [**canopy\_date\_mod::nfeb**](#function-nfeb) (year year) <br> |
|  subroutine | [**canopy\_date\_mod::split\_date\_char**](#function-split_date_char) (date date, century\_year century\_year, month month, day day, hour hour, minute minute, second second) <br> |




























## Detailed Description


This group contains subroutines and functions for:
* Computing time differences between dates
* Generating new dates from time increments
* Converting between date formats
* Julian day calculations
* Date string parsing and manipulation 




    
## Public Functions Documentation




### function geth\_idts 

_Calculate time difference between two dates._ 
```Fortran
subroutine canopy_date_mod::geth_idts (
    newdate newdate,
    olddate olddate,
    idts idts
) 
```



Computes the time difference in seconds between two input dates in the format 'YYYY-MM-DD HH:MM:SS.ffff'. Handles leap years and different month lengths correctly. 

**Parameters:**


* `newdate` The new (later) date string 'YYYY-MM-DD HH:MM:SS.ffff' 
* `olddate` The old (earlier) date string 'YYYY-MM-DD HH:MM:SS.ffff' 
* `idts` Time difference in seconds (newdate - olddate) [output] 



**Author:**

NCAR, T. Otte 




**Date:**

February 2001, modified 2011 





        

<hr>



### function geth\_newdate 

```Fortran
subroutine canopy_date_mod::geth_newdate (
    ndate ndate,
    odate odate,
    idt idt
) 
```




<hr>



### function getsdt 

```Fortran
subroutine canopy_date_mod::getsdt (
    hdate hdate,
    sdate sdate,
    stime stime
) 
```




<hr>



### function julian 

```Fortran
integer  function canopy_date_mod::julian (
    year year,
    mnth mnth,
    mday mday
) 
```




<hr>



### function nfeb 

```Fortran
integer  function canopy_date_mod::nfeb (
    year year
) 
```




<hr>



### function split\_date\_char 

```Fortran
subroutine canopy_date_mod::split_date_char (
    date date,
    century_year century_year,
    month month,
    day day,
    hour hour,
    minute minute,
    second second
) 
```




<hr>

------------------------------


