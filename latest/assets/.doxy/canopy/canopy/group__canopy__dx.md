

# Group canopy\_dx



[**Modules**](modules.md) **>** [**canopy\_dx**](group__canopy__dx.md)



_Grid cell distance calculations using Haversine formula._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_dxcalc\_mod::canopy\_calcdx**](#function-canopy_calcdx) (DXOPT DXOPT, DXSET DXSET, NLAT NLAT, NLON NLON, LAT LAT, LON LON, DX DX) <br>_Calculate grid cell distances for 1D arrays._  |
|  subroutine | [**canopy\_dxcalc\_mod::canopy\_calcdx\_2d**](#function-canopy_calcdx_2d) (DXOPT DXOPT, DXSET DXSET, NLAT NLAT, NLON NLON, LAT LAT, LON LON, DX DX) <br>_Calculate grid cell distances for 2D arrays._  |




























## Public Functions Documentation




### function canopy\_calcdx 

_Calculate grid cell distances for 1D arrays._ 
```Fortran
subroutine canopy_dxcalc_mod::canopy_calcdx (
    DXOPT DXOPT,
    DXSET DXSET,
    NLAT NLAT,
    NLON NLON,
    LAT LAT,
    LON LON,
    DX DX
) 
```



Computes great circle distance or the orthodromic distance using Haversine formula for 1D latitude/longitude arrays




**Parameters:**


* `DXOPT` User DX calculation option (0=calculate, 1=use set value) 
* `DXSET` User DX set value if cannot calculate (m) 
* `NLAT` Number of latitude grid cells/points 
* `NLON` Number of longitude grid cells/points 
* `LAT` Model latitudes array (degrees) 
* `LON` Model longitudes array (degrees) 
* `DX` Distance between two points (m) 



**Author:**

P. C. Campbell 




**Date:**

Oct 2022 




**Note:**

Uses Haversine formula for great circle distance calculations 





        

<hr>



### function canopy\_calcdx\_2d 

_Calculate grid cell distances for 2D arrays._ 
```Fortran
subroutine canopy_dxcalc_mod::canopy_calcdx_2d (
    DXOPT DXOPT,
    DXSET DXSET,
    NLAT NLAT,
    NLON NLON,
    LAT LAT,
    LON LON,
    DX DX
) 
```



Computes great circle distance or the orthodromic distance using Haversine formula for 2D latitude/longitude arrays




**Parameters:**


* `DXOPT` User DX calculation option (0=calculate, 1=use set value) 
* `DXSET` User DX set value if cannot calculate (m) 
* `NLAT` Number of latitude grid cells/points 
* `NLON` Number of longitude grid cells/points 
* `LAT` Model latitudes 2D array (degrees) 
* `LON` Model longitudes 2D array (degrees) 
* `DX` Distance between two points 2D array (m) 



**Author:**

P. C. Campbell 




**Date:**

Oct 2022 




**Note:**

Uses Haversine formula for great circle distance calculations 





        

<hr>

------------------------------


