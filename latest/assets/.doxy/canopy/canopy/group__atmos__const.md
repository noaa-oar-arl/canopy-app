

# Group atmos\_const



[**Modules**](modules.md) **>** [**atmos\_const**](group__atmos__const.md)



_Constants for atmospheric thermodynamics and chemistry._ 






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::beta\_n**](#variable-beta_n)   = `0.35\_rk`<br>_Stability parameter beta for neutral conditions (dimensionless)_  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::cpd**](#variable-cpd)   = `7.0 \* [**rdgas**](namespacecanopy__const__mod.md#variable-rdgas) / 2.0`<br>_Specific heat of dry air at constant pressure [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::cpwvap**](#variable-cpwvap)   = `4.0 \* [**rwvap**](namespacecanopy__const__mod.md#variable-rwvap)`<br>_Specific heat for water vapor at constant pressure [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::cvd**](#variable-cvd)   = `5.0 \* [**rdgas**](namespacecanopy__const__mod.md#variable-rdgas) / 2.0`<br>_Specific heat of dry air at constant volume [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::cvwvap**](#variable-cvwvap)   = `3.0 \* [**rwvap**](namespacecanopy__const__mod.md#variable-rwvap)`<br>_Specific heat for water vapor at constant volume [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::dlvdt**](#variable-dlvdt)   = `2370.0`<br>_Rate of change of latent heat of vaporization w.r.t._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::lf0**](#variable-lf0)   = `3.34e5`<br>_Latent heat of fusion of water at 0°C [J/kg]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::lv0**](#variable-lv0)   = `2.501e6`<br>_Latent heat of vaporization of water at 0°C [J/kg]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::mwair**](#variable-mwair)   = `28.9628\_rk`<br>_Mean molecular weight for dry air [g/mol]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::mwwat**](#variable-mwwat)   = `18.01528\_rk`<br>_Mean molecular weight for water vapor [g/mol]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::rdgas**](#variable-rdgas)   = `1.0e3 \* [**rgasuniv**](namespacecanopy__const__mod.md#variable-rgasuniv) / [**mwair**](namespacecanopy__const__mod.md#variable-mwair)`<br>_Dry-air gas constant [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::rwvap**](#variable-rwvap)   = `1.0e3 \* [**rgasuniv**](namespacecanopy__const__mod.md#variable-rgasuniv) / [**mwwat**](namespacecanopy__const__mod.md#variable-mwwat)`<br>_Gas constant for water vapor [J/kg-K]._  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::vonk**](#variable-vonk)   = `0.4\_rk`<br>_Von Karman constant (dimensionless)_  |
|  real([**rk**](namespacecanopy__const__mod.md#variable-rk)), parameter | [**canopy\_const\_mod::vp0**](#variable-vp0)   = `611.29\_rk`<br>_Vapor pressure of water at 0°C_ [Pa](Source: CRC76 pp. 6-15) __ |












































## Public Attributes Documentation




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



### variable dlvdt 

_Rate of change of latent heat of vaporization w.r.t._ 
```Fortran
real(rk), parameter canopy_const_mod::dlvdt;
```



temperature [J/kg-K] 


        

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



### variable rdgas 

_Dry-air gas constant [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::rdgas;
```



287.07548994 J/kg-K 


        

<hr>



### variable rwvap 

_Gas constant for water vapor [J/kg-K]._ 
```Fortran
real(rk), parameter canopy_const_mod::rwvap;
```



461.52492604 J/kg-K 


        

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


