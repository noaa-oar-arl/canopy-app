

# Namespace canopy\_const\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_const\_mod**](namespacecanopy__const__mod.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**avo**](#variable-avo)   = `6.02214076d23`<br>_Avogadro's Constant [number/mol]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**beta\_n**](#variable-beta_n)   = `0.35\_rk`<br>_Stability parameter beta for neutral conditions (dimensionless)_  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**cpd**](#variable-cpd)   = `7.0 \* [**rdgas**](namespacecanopy__const__mod.md#variable-rdgas) / 2.0`<br>_Specific heat of dry air at constant pressure [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**cpwvap**](#variable-cpwvap)   = `4.0 \* [**rwvap**](namespacecanopy__const__mod.md#variable-rwvap)`<br>_Specific heat for water vapor at constant pressure [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**cvd**](#variable-cvd)   = `5.0 \* [**rdgas**](namespacecanopy__const__mod.md#variable-rdgas) / 2.0`<br>_Specific heat of dry air at constant volume [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**cvwvap**](#variable-cvwvap)   = `3.0 \* [**rwvap**](namespacecanopy__const__mod.md#variable-rwvap)`<br>_Specific heat for water vapor at constant volume [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)) | [**dg2m**](#variable-dg2m)  <br>_Latitude degrees to meters conversion factor [m/deg]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**dlvdt**](#variable-dlvdt)   = `2370.0`<br>_Rate of change of latent heat of vaporization w.r.t._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**fillreal**](#variable-fillreal)   = `-9.0e20`<br>_NetCDF fill value for missing real data._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**grav**](#variable-grav)   = `9.80622\_rk`<br>_Mean gravitational acceleration [m/sec²]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**lf0**](#variable-lf0)   = `3.34e5`<br>_Latent heat of fusion of water at 0°C [J/kg]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**lv0**](#variable-lv0)   = `2.501e6`<br>_Latent heat of vaporization of water at 0°C [J/kg]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**molvol**](#variable-molvol)   = `22.4139695\_rk`<br>_Molar volume of ideal gas at STP [L/mol] (Non-MKS units)_  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**mwair**](#variable-mwair)   = `28.9628\_rk`<br>_Mean molecular weight for dry air [g/mol]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**mwwat**](#variable-mwwat)   = `18.01528\_rk`<br>_Mean molecular weight for water vapor [g/mol]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**pi**](#variable-pi)   = `3.14159265358979324\_rk`<br>_Pi (double precision: 3.14159265358979324)_  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**pi180**](#variable-pi180)   = `[**pi**](namespacecanopy__const__mod.md#variable-pi) / 180.0\_rk`<br>_Pi/180 conversion factor from degrees to radians [rad/deg]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**rdgas**](#variable-rdgas)   = `1.0e3 \* [**rgasuniv**](namespacecanopy__const__mod.md#variable-rgasuniv) / [**mwair**](namespacecanopy__const__mod.md#variable-mwair)`<br>_Dry-air gas constant [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**rearth**](#variable-rearth)   = `6371008.8\_rk`<br>_Radius of Earth [m]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**rgasuniv**](#variable-rgasuniv)   = `8.31446261815324\_rk`<br>_Universal gas constant [J/mol-K]._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), parameter | [**rk**](#variable-rk)   = `SELECTED\_REAL\_KIND(15, 307)`<br>_Selected real kind for double precision (15 decimal digits, 307 exponent range)_  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**rwvap**](#variable-rwvap)   = `1.0e3 \* [**rgasuniv**](namespacecanopy__const__mod.md#variable-rgasuniv) / [**mwwat**](namespacecanopy__const__mod.md#variable-mwwat)`<br>_Gas constant for water vapor [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**siday**](#variable-siday)   = `86164.09\_rk`<br>_Length of a sidereal day_ [sec](Source: CRC76, pp. 14-6) __ |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**solcnst**](#variable-solcnst)   = `1373.0`<br>_Solar Constant_ [W/m²](Source: CRC76, pp. 14-2) __ |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**stdatmpa**](#variable-stdatmpa)   = `101325.0`<br>_Standard atmosphere [Pa]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**stdtemp**](#variable-stdtemp)   = `273.15\_rk`<br>_Standard temperature [K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**stfblz**](#variable-stfblz)   = `5.67037442[**d**](canopy__calcs_8F90.md#variable-d)-8`<br>_Stefan-Boltzmann constant [W/(m² K⁴)]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**tau\_days**](#variable-tau_days)   = `5.0\_rk`<br>_E-folding time for long-term past conditions for biogenic emissions [days]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**tau\_hours**](#variable-tau_hours)   = `12.0\_rk`<br>_E-folding time for short-term past conditions for biogenic emissions [hours]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**vonk**](#variable-vonk)   = `0.4\_rk`<br>_Von Karman constant (dimensionless)_  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**vp0**](#variable-vp0)   = `611.29\_rk`<br>_Vapor pressure of water at 0°C_ [Pa](Source: CRC76 pp. 6-15) __ |












































## Public Attributes Documentation




### variable avo 

_Avogadro's Constant [number/mol]._ 
```Fortran
real(rk), parameter canopy_const_mod::avo;
```




<hr>



### variable beta\_n 

_Stability parameter beta for neutral conditions (dimensionless)_ 
```Fortran
real(rk), parameter canopy_const_mod::beta_n;
```



Source: Bonan et al. (2018) [https://doi.org/10.5194/gmd-11-1467-2018](https://doi.org/10.5194/gmd-11-1467-2018) 


        

<hr>



### variable cpd 

_Specific heat of dry air at constant pressure [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::cpd;
```



1004.7642148 J/kg-K. Calculated assuming dry air is classical ideal gas 


        

<hr>



### variable cpwvap 

_Specific heat for water vapor at constant pressure [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::cpwvap;
```



1846.0997042 J/kg-K. Calculated assuming water vapor is classical ideal gas 


        

<hr>



### variable cvd 

_Specific heat of dry air at constant volume [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::cvd;
```



717.68872485 J/kg-K. Calculated assuming dry air is classical ideal gas 


        

<hr>



### variable cvwvap 

_Specific heat for water vapor at constant volume [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::cvwvap;
```



1384.5747781 J/kg-K. Calculated assuming water vapor is classical ideal gas 


        

<hr>



### variable dg2m 

_Latitude degrees to meters conversion factor [m/deg]._ 
```Fortran
real(rk) canopy_const_mod::dg2m;
```




<hr>



### variable dlvdt 

_Rate of change of latent heat of vaporization w.r.t._ 
```Fortran
real(rk), parameter canopy_const_mod::dlvdt;
```



temperature [J/kg-K] 


        

<hr>



### variable fillreal 

_NetCDF fill value for missing real data._ 
```Fortran
real(rk), parameter canopy_const_mod::fillreal;
```




<hr>



### variable grav 

_Mean gravitational acceleration [m/sec²]._ 
```Fortran
real(rk), parameter canopy_const_mod::grav;
```



Mean of polar and equatorial values (Source: CRC76, pp. 14-6) 


        

<hr>



### variable lf0 

_Latent heat of fusion of water at 0°C [J/kg]._ 
```Fortran
real(rk), parameter canopy_const_mod::lf0;
```




<hr>



### variable lv0 

_Latent heat of vaporization of water at 0°C [J/kg]._ 
```Fortran
real(rk), parameter canopy_const_mod::lv0;
```



Values from p. 641 of Stull (1988) 


        

<hr>



### variable molvol 

_Molar volume of ideal gas at STP [L/mol] (Non-MKS units)_ 
```Fortran
real(rk), parameter canopy_const_mod::molvol;
```




<hr>



### variable mwair 

_Mean molecular weight for dry air [g/mol]._ 
```Fortran
real(rk), parameter canopy_const_mod::mwair;
```



78.06% N2, 21% O2, and 0.943% A on a mole fraction basis (Source: Hobbs, 1995, pp. 69-70) 


        

<hr>



### variable mwwat 

_Mean molecular weight for water vapor [g/mol]._ 
```Fortran
real(rk), parameter canopy_const_mod::mwwat;
```




<hr>



### variable pi 

_Pi (double precision: 3.14159265358979324)_ 
```Fortran
real(rk), parameter canopy_const_mod::pi;
```




<hr>



### variable pi180 

_Pi/180 conversion factor from degrees to radians [rad/deg]._ 
```Fortran
real(rk), parameter canopy_const_mod::pi180;
```




<hr>



### variable rdgas 

_Dry-air gas constant [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::rdgas;
```



287.07548994 J/kg-K 


        

<hr>



### variable rearth 

_Radius of Earth [m]._ 
```Fortran
real(rk), parameter canopy_const_mod::rearth;
```



Radius of sphere having same surface area as Clarke ellipsoid of 1866 WGS84 arithmetic mean radius (Source: Snyder, 1987) 


        

<hr>



### variable rgasuniv 

_Universal gas constant [J/mol-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::rgasuniv;
```




<hr>



### variable rk 

_Selected real kind for double precision (15 decimal digits, 307 exponent range)_ 
```Fortran
integer, parameter canopy_const_mod::rk;
```




<hr>



### variable rwvap 

_Gas constant for water vapor [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::rwvap;
```



461.52492604 J/kg-K 


        

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



### variable stdatmpa 

_Standard atmosphere [Pa]._ 
```Fortran
real(rk), parameter canopy_const_mod::stdatmpa;
```




<hr>



### variable stdtemp 

_Standard temperature [K]._ 
```Fortran
real(rk), parameter canopy_const_mod::stdtemp;
```




<hr>



### variable stfblz 

_Stefan-Boltzmann constant [W/(m² K⁴)]._ 
```Fortran
real(rk), parameter canopy_const_mod::stfblz;
```




<hr>



### variable tau\_days 

_E-folding time for long-term past conditions for biogenic emissions [days]._ 
```Fortran
real(rk), parameter canopy_const_mod::tau_days;
```




<hr>



### variable tau\_hours 

_E-folding time for short-term past conditions for biogenic emissions [hours]._ 
```Fortran
real(rk), parameter canopy_const_mod::tau_hours;
```




<hr>



### variable vonk 

_Von Karman constant (dimensionless)_ 
```Fortran
real(rk), parameter canopy_const_mod::vonk;
```




<hr>



### variable vp0 

_Vapor pressure of water at 0°C_ [Pa](Source: CRC76 pp. 6-15) __
```Fortran
real(rk), parameter canopy_const_mod::vp0;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_const_mod.F90`

