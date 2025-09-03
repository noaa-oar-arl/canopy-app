

# Group geodetic\_const



[**Modules**](modules.md) **>** [**geodetic\_const**](group__geodetic__const.md)



_Earth geometry and geodetic constants._ 






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)) | [**canopy\_const\_mod::dg2m**](#variable-dg2m)  <br>_Latitude degrees to meters conversion factor [m/deg]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::grav**](#variable-grav)   = `9.80622\_rk`<br>_Mean gravitational acceleration [m/sec²]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::rearth**](#variable-rearth)   = `6371008.8\_rk`<br>_Radius of Earth [m]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::siday**](#variable-siday)   = `86164.09\_rk`<br>_Length of a sidereal day_ [sec](Source: CRC76, pp. 14-6) __ |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::solcnst**](#variable-solcnst)   = `1373.0`<br>_Solar Constant_ [W/m²](Source: CRC76, pp. 14-2) __ |












































## Public Attributes Documentation




### variable dg2m 

_Latitude degrees to meters conversion factor [m/deg]._ 
```Fortran
real(rk) canopy_const_mod::dg2m;
```




<hr>



### variable grav 

_Mean gravitational acceleration [m/sec²]._ 
```Fortran
real(rk), parameter canopy_const_mod::grav;
```



Mean of polar and equatorial values (Source: CRC76, pp. 14-6) 


        

<hr>



### variable rearth 

_Radius of Earth [m]._ 
```Fortran
real(rk), parameter canopy_const_mod::rearth;
```



Radius of sphere having same surface area as Clarke ellipsoid of 1866 WGS84 arithmetic mean radius (Source: Snyder, 1987) 


        

<hr>



### variable siday 

_Length of a sidereal day_ [sec](Source: CRC76, pp. 14-6) __
```Fortran
real(rk), parameter canopy_const_mod::siday;
```




<hr>



### variable solcnst 

_Solar Constant_ [W/m²](Source: CRC76, pp. 14-2) __
```Fortran
real(rk), parameter canopy_const_mod::solcnst;
```




<hr>

------------------------------


