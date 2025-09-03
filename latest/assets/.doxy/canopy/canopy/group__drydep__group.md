

# Group drydep\_group



[**Modules**](modules.md) **>** [**drydep\_group**](group__drydep__group.md)



_Routines for calculating gas dry deposition to various surface types._ [More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  real(rk), parameter | [**ar\_0**](#variable-ar_0)   = `8.0`<br>_Used to scale other species to HNO3 (dimensionless)_  |
|  real(rk) | [**ar\_l**](#variable-ar_l)  <br>_Reactivity denominator relative to HNO3 for each species (dimensionless)_  |
|  real(rk) | [**cave\_l**](#variable-cave_l)  <br>_Maxwell-Boltzmann average speed of gas distribution (m/s)_  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**chemmechgas\_opt**](#variable-chemmechgas_opt)  <br>_Select chemical mechanism._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**chemmechgas\_tot**](#variable-chemmechgas_tot)  <br>_Select chemical mechanism gas species list._  |
|  real(rk) | [**cp\_air**](#variable-cp_air)  <br>_Specific heat of moist air (J/kg-K)_  |
|  real(rk) | [**ctemp2**](#variable-ctemp2)  <br>_Mean temperature just above surface (C)_  |
|  real(rk), parameter | [**d3**](#variable-d3)   = `1.38564e-2`<br>_Scaling parameter used to estimate the friction velocity in surface waters from the atmospheric friction velocity to a value following Slinn et al._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**dep\_ind**](#variable-dep_ind)  <br>_Gas deposition species index (depends on gas mech, set in constants)_  |
|  real(rk), dimension(:), intent([**out**](canopy__bioparm__mod_8F90.md#variable-out)) | [**dep\_out**](#variable-dep_out)  <br>_Output canopy layer gas dry deposition rate for each DEP\_IND (cm/s)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**dsoil**](#variable-dsoil)  <br>_Depth of topsoil (cm)_  |
|  real(rk) | [**dw**](#variable-dw)  <br>_Diffusivity of water._  |
|  real(rk) | [**dw25**](#variable-dw25)  <br>_Diffusivity of water at 298.15 K._  |
|  real(rk) | [**f01**](#variable-f01)  <br>_Reactivity parameter based on DEP\_IND (0-1)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fch**](#variable-fch)  <br>_Canopy height (m)_  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**fsun**](#variable-fsun)  <br>_Sunlit/Shaded fraction from photolysis correction factor._  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**gamma\_build**](#variable-gamma_build)  <br>_Reaction probability with building type (dimensionless) Default NL is average of range in gamma from as low as 10−8 for glass and metal to 10−4 for activated carbon and brick._  |
|  real(rk) | [**hstarl**](#variable-hstarl)  <br>_Effective Henry's law coefficient based on DEP\_IND (M/atm)_  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer) | [**i**](#variable-i)  <br>_Loop index._  |
|  real(rk) | [**kvisw**](#variable-kvisw)  <br>_Kinematic viscosity of water (cm^2/s)_  |
|  real(rk) | [**lebas\_l**](#variable-lebas_l)  <br>_Le Bas molar volumes are from the Schroeder additive method (cm3/mol)_  |
|  real(rk) | [**lv**](#variable-lv)  <br>_Latent heat of vaporization (J/kg)_  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**mdiffl**](#variable-mdiffl)  <br>_Molecular diffusivity for species l based on DEP\_IND (cm2/s)_  |
|  real(rk) | [**mmg\_l**](#variable-mmg_l)  <br>_Molar mass for each gas species (kg/mol)_  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**ppfd**](#variable-ppfd)  <br>_PPFD ave sun and shade (umol/m2 s)_  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd\_shade**](#variable-ppfd_shade)  <br>_PPFD for shaded leaves (umol phot/m2 s)_  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ppfd\_sun**](#variable-ppfd_sun)  <br>_PPFD for sunlit leaves (umol phot/m2 s)_  |
|  real(rk), parameter | [**pr**](#variable-pr)   = `0.709`<br>_Prandtl Number (dimensionless)_  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**pressa**](#variable-pressa)  <br>_Ambient Pressure profile in canopy (mb)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**qv2**](#variable-qv2)  <br>_Mean mixing ratio just above surface (kg/kg)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ra**](#variable-ra)  <br>_Aerodynamic resistance (s/cm)_  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**rb**](#variable-rb)  <br>_Boundary layer resistance for species l based on DEP\_IND (s/cm)_  |
|  real(rk) | [**rbg**](#variable-rbg)  <br>_Ground boundary layer resistance (s/cm)_  |
|  real(rk) | [**rbw**](#variable-rbw)  <br>_Water boundary layer resistance (s/cm)_  |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**rc**](#variable-rc)  <br>_Cuticular resistance for species l based on DEP\_IND (s/cm)_  |
|  real(rk) | [**rden**](#variable-rden)  <br> |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**relhuma**](#variable-relhuma)  <br>_Ambient Relative Humidity profile in canopy (%)_  |
|  real(rk) | [**rlx**](#variable-rlx)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**rm**](#variable-rm)  <br>_Mesophyll resistance for species l based on DEP\_IND (s/cm)_  |
|  real(rk) | [**rnum**](#variable-rnum)  <br> |
|  real(rk), dimension(size([**zk**](canopy__var3din__mod_8F90.md#variable-zk))) | [**rs**](#variable-rs)  <br>_Stomatal resistance for species l based on DEP\_IND (s/cm)_  |
|  real(rk), parameter | [**rsnow0**](#variable-rsnow0)   = `100.0`<br>_Resistance to deposition to snow (s/cm) based on Helmig et al._  |
|  real(rk) | [**rsnowl**](#variable-rsnowl)  <br>_Resistance to diffusion thru snow space for chemical species (s/cm)_  |
|  real(rk) | [**rsoill**](#variable-rsoill)  <br>_Resistance to diffusion thru soil pore space for chemical species (s/cm)_  |
|  real(rk), parameter | [**rt25ink**](#variable-rt25ink)   = `1.0\_rk/(stdtemp + 25.0\_rk)`<br>_298.15K = 25C_  |
|  real(rk) | [**rurbanl**](#variable-rurbanl)  <br>_Resistance to diffusion thru snow space for chemical species (s/m)_  |
|  real(rk) | [**rwaterl**](#variable-rwaterl)  <br>_Water surface resistance (s/cm)_  |
|  real(rk) | [**scw\_pr\_23**](#variable-scw_pr_23)  <br>_(scw/pr)\*\*2/3_  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**socat**](#variable-socat)  <br>_Input soil category dataset used._  |
|  [**integer**](canopy__bioparm__mod_8F90.md#variable-integer), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**sotyp**](#variable-sotyp)  <br>_Input soil type integer associated with soilcat._  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**srad**](#variable-srad)  <br>_Incoming solar irradiation top of canopy (W/m^2)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**stheta**](#variable-stheta)  <br>_Volumetric soil water content in topsoil (m^3/m^3)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**temp**](#variable-temp)  <br>_Mean temperature just above surface (K)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**temp2**](#variable-temp2)  <br>_Mean temperature just above surface (K)_  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tempa**](#variable-tempa)  <br>_Ambient Temperature profile in canopy (K)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**tempsoil**](#variable-tempsoil)  <br>_Soil temperature in topsoil (K)_  |
|  real(rk) | [**tw**](#variable-tw)  <br>_Wet bulb temperature (K)_  |
|  real(rk), parameter | [**twothirds**](#variable-twothirds)   = `2.0\_rk / 3.0\_rk`<br>_Two thirds constant._  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ubar**](#variable-ubar)  <br>_Mean above/in-canopy wind speed (m/s)_  |
|  real(rk), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**ustar**](#variable-ustar)  <br>_Friction velocity at surface (m/s)_  |
|  real(rk) | [**vdlx**](#variable-vdlx)  <br>_Working variables for resistance calculations._  |
|  real(rk), dimension(:), intent([**in**](canopy__calcs_8F90.md#variable-in)) | [**zk**](#variable-zk)  <br>_Model heights (m)_  |












































## Detailed Description


This group contains subroutines for calculating gas dry deposition velocities to vegetation, soil, snow, urban surfaces, and water based on Zhang et al. (2003) and other established parameterizations. 


    
## Public Attributes Documentation




### variable ar\_0 

_Used to scale other species to HNO3 (dimensionless)_ 
```Fortran
real(rk), parameter ar_0;
```




<hr>



### variable ar\_l 

_Reactivity denominator relative to HNO3 for each species (dimensionless)_ 
```Fortran
real(rk) ar_l;
```




<hr>



### variable cave\_l 

_Maxwell-Boltzmann average speed of gas distribution (m/s)_ 
```Fortran
real(rk) cave_l;
```




<hr>



### variable chemmechgas\_opt 

_Select chemical mechanism._ 
```Fortran
integer, intent(in) chemmechgas_opt;
```




<hr>



### variable chemmechgas\_tot 

_Select chemical mechanism gas species list._ 
```Fortran
integer, intent(in) chemmechgas_tot;
```




<hr>



### variable cp\_air 

_Specific heat of moist air (J/kg-K)_ 
```Fortran
real(rk) cp_air;
```




<hr>



### variable ctemp2 

_Mean temperature just above surface (C)_ 
```Fortran
real(rk) ctemp2;
```




<hr>



### variable d3 

_Scaling parameter used to estimate the friction velocity in surface waters from the atmospheric friction velocity to a value following Slinn et al._ 
```Fortran
real(rk), parameter d3;
```



(1978) and Fairall et al. (2007) 


        

<hr>



### variable dep\_ind 

_Gas deposition species index (depends on gas mech, set in constants)_ 
```Fortran
integer, intent(in) dep_ind;
```



Gas deposition species index (depends on gas mech) 


        

<hr>



### variable dep\_out 

_Output canopy layer gas dry deposition rate for each DEP\_IND (cm/s)_ 
```Fortran
real(rk), intent(out) dep_out;
```



Output soil layer gas dry deposition rate for each DEP\_IND (cm/s) 


        

<hr>



### variable dsoil 

_Depth of topsoil (cm)_ 
```Fortran
real(rk), intent(in) dsoil;
```




<hr>



### variable dw 

_Diffusivity of water._ 
```Fortran
real(rk) dw;
```




<hr>



### variable dw25 

_Diffusivity of water at 298.15 K._ 
```Fortran
real(rk) dw25;
```




<hr>



### variable f01 

_Reactivity parameter based on DEP\_IND (0-1)_ 
```Fortran
real(rk) f01;
```




<hr>



### variable fch 

_Canopy height (m)_ 
```Fortran
real(rk), intent(in) fch;
```




<hr>



### variable fsun 

_Sunlit/Shaded fraction from photolysis correction factor._ 
```Fortran
real(rk), dimension(:), intent(in) fsun;
```




<hr>



### variable gamma\_build 

_Reaction probability with building type (dimensionless) Default NL is average of range in gamma from as low as 10−8 for glass and metal to 10−4 for activated carbon and brick._ 
```Fortran
real(rk), intent(in) gamma_build;
```



=5.0D-5. Reference (Shen and Gao, 2018; [https://doi.org/10.1016/j.buildenv.2018.02.046](https://doi.org/10.1016/j.buildenv.2018.02.046)) 


        

<hr>



### variable hstarl 

_Effective Henry's law coefficient based on DEP\_IND (M/atm)_ 
```Fortran
real(rk) hstarl;
```




<hr>



### variable i 

_Loop index._ 
```Fortran
integer i;
```




<hr>



### variable kvisw 

_Kinematic viscosity of water (cm^2/s)_ 
```Fortran
real(rk) kvisw;
```




<hr>



### variable lebas\_l 

_Le Bas molar volumes are from the Schroeder additive method (cm3/mol)_ 
```Fortran
real(rk) lebas_l;
```




<hr>



### variable lv 

_Latent heat of vaporization (J/kg)_ 
```Fortran
real(rk) lv;
```




<hr>



### variable mdiffl 

_Molecular diffusivity for species l based on DEP\_IND (cm2/s)_ 
```Fortran
real(rk) mdiffl;
```



Molecular diffusivity (cm^2/s) 


        

<hr>



### variable mmg\_l 

_Molar mass for each gas species (kg/mol)_ 
```Fortran
real(rk) mmg_l;
```




<hr>



### variable ppfd 

_PPFD ave sun and shade (umol/m2 s)_ 
```Fortran
real(rk), dimension(size(zk)) ppfd;
```




<hr>



### variable ppfd\_shade 

_PPFD for shaded leaves (umol phot/m2 s)_ 
```Fortran
real(rk), dimension(:), intent(in) ppfd_shade;
```




<hr>



### variable ppfd\_sun 

_PPFD for sunlit leaves (umol phot/m2 s)_ 
```Fortran
real(rk), dimension(:), intent(in) ppfd_sun;
```




<hr>



### variable pr 

_Prandtl Number (dimensionless)_ 
```Fortran
real(rk), parameter pr;
```




<hr>



### variable pressa 

_Ambient Pressure profile in canopy (mb)_ 
```Fortran
real(rk), intent(in) pressa;
```



Ambient Pressure just above surface (mb) 


        

<hr>



### variable qv2 

_Mean mixing ratio just above surface (kg/kg)_ 
```Fortran
real(rk), intent(in) qv2;
```




<hr>



### variable ra 

_Aerodynamic resistance (s/cm)_ 
```Fortran
real(rk), intent(in) ra;
```




<hr>



### variable rb 

_Boundary layer resistance for species l based on DEP\_IND (s/cm)_ 
```Fortran
real(rk), dimension(size(zk)) rb;
```




<hr>



### variable rbg 

_Ground boundary layer resistance (s/cm)_ 
```Fortran
real(rk) rbg;
```




<hr>



### variable rbw 

_Water boundary layer resistance (s/cm)_ 
```Fortran
real(rk) rbw;
```




<hr>



### variable rc 

_Cuticular resistance for species l based on DEP\_IND (s/cm)_ 
```Fortran
real(rk), dimension(size(zk)) rc;
```




<hr>



### variable rden 

```Fortran
real(rk) rden;
```




<hr>



### variable relhuma 

_Ambient Relative Humidity profile in canopy (%)_ 
```Fortran
real(rk), dimension(:), intent(in) relhuma;
```




<hr>



### variable rlx 

```Fortran
real(rk) rlx;
```




<hr>



### variable rm 

_Mesophyll resistance for species l based on DEP\_IND (s/cm)_ 
```Fortran
real(rk), dimension(size(zk)) rm;
```




<hr>



### variable rnum 

```Fortran
real(rk) rnum;
```




<hr>



### variable rs 

_Stomatal resistance for species l based on DEP\_IND (s/cm)_ 
```Fortran
real(rk), dimension(size(zk)) rs;
```




<hr>



### variable rsnow0 

_Resistance to deposition to snow (s/cm) based on Helmig et al._ 
```Fortran
real(rk), parameter rsnow0;
```




<hr>



### variable rsnowl 

_Resistance to diffusion thru snow space for chemical species (s/cm)_ 
```Fortran
real(rk) rsnowl;
```




<hr>



### variable rsoill 

_Resistance to diffusion thru soil pore space for chemical species (s/cm)_ 
```Fortran
real(rk) rsoill;
```




<hr>



### variable rt25ink 

_298.15K = 25C_ 
```Fortran
real(rk), parameter rt25ink;
```




<hr>



### variable rurbanl 

_Resistance to diffusion thru snow space for chemical species (s/m)_ 
```Fortran
real(rk) rurbanl;
```




<hr>



### variable rwaterl 

_Water surface resistance (s/cm)_ 
```Fortran
real(rk) rwaterl;
```




<hr>



### variable scw\_pr\_23 

_(scw/pr)\*\*2/3_ 
```Fortran
real(rk) scw_pr_23;
```




<hr>



### variable socat 

_Input soil category dataset used._ 
```Fortran
integer, intent(in) socat;
```




<hr>



### variable sotyp 

_Input soil type integer associated with soilcat._ 
```Fortran
integer, intent(in) sotyp;
```




<hr>



### variable srad 

_Incoming solar irradiation top of canopy (W/m^2)_ 
```Fortran
real(rk), intent(in) srad;
```




<hr>



### variable stheta 

_Volumetric soil water content in topsoil (m^3/m^3)_ 
```Fortran
real(rk), intent(in) stheta;
```




<hr>



### variable temp 

_Mean temperature just above surface (K)_ 
```Fortran
real(rk), intent(in) temp;
```




<hr>



### variable temp2 

_Mean temperature just above surface (K)_ 
```Fortran
real(rk), intent(in) temp2;
```




<hr>



### variable tempa 

_Ambient Temperature profile in canopy (K)_ 
```Fortran
real(rk), dimension(:), intent(in) tempa;
```




<hr>



### variable tempsoil 

_Soil temperature in topsoil (K)_ 
```Fortran
real(rk), intent(in) tempsoil;
```




<hr>



### variable tw 

_Wet bulb temperature (K)_ 
```Fortran
real(rk) tw;
```




<hr>



### variable twothirds 

_Two thirds constant._ 
```Fortran
real(rk), parameter twothirds;
```




<hr>



### variable ubar 

_Mean above/in-canopy wind speed (m/s)_ 
```Fortran
real(rk), intent(in) ubar;
```



Mean wind speed just above surface (m/s) 


        

<hr>



### variable ustar 

_Friction velocity at surface (m/s)_ 
```Fortran
real(rk), intent(in) ustar;
```




<hr>



### variable vdlx 

_Working variables for resistance calculations._ 
```Fortran
real(rk) vdlx;
```




<hr>



### variable zk 

_Model heights (m)_ 
```Fortran
real(rk), dimension(:), intent(in) zk;
```




<hr>

------------------------------


