

# Namespace canopy\_aero\_ddep\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_aero\_ddep\_mod**](namespacecanopy__aero__ddep__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine, public | [**canopy\_aero\_ddep\_pleim2022**](#function-canopy_aero_ddep_pleim2022) ([**ustar**](canopy__drydep__mod_8F90.md#variable-ustar) ustar, Ra Ra, u u, d\_p d\_p, rho\_p rho\_p, T T, P P, surface\_type surface\_type, vdep vdep) <br>_Aerosol dry deposition velocity for non-vegetated/urban areas using Pleim et al._  |
|  subroutine, public | [**canopy\_aero\_ddep\_subveg**](#function-canopy_aero_ddep_subveg) (nlev nlev, [**z**](canopy__var3din__mod_8F90.md#variable-z) z, hc hc, lad lad, u u, d\_p d\_p, rho\_p rho\_p, T T, P P, vdep\_opt vdep\_opt, Ra Ra, [**modres**](canopy__bioemi__mod_8F90.md#variable-modres) modres, [**ustar**](canopy__drydep__mod_8F90.md#variable-ustar) ustar, [**vtype**](canopy__profile__mod_8F90.md#variable-vtype) vtype, vdep vdep) <br>_Sub-canopy aerosol dry deposition velocity following Pleim et al._  |




























## Public Functions Documentation




### function canopy\_aero\_ddep\_pleim2022 

_Aerosol dry deposition velocity for non-vegetated/urban areas using Pleim et al._ 
```Fortran
subroutine, public canopy_aero_ddep_mod::canopy_aero_ddep_pleim2022 (
    ustar ustar,
    Ra Ra,
    u u,
    d_p d_p,
    rho_p rho_p,
    T T,
    P P,
    surface_type surface_type,
    vdep vdep
) 
```



(2022)


Implements Pleim et al. (2022) equations for urban and non-vegetated surfaces. [https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022MS003050](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022MS003050) 

**Parameters:**


* `u` Array of wind speed (m/s) 
* `d_p` Aerosol particle diameter (m) 
* `rho_p` Particle density (kg/m^3) 
* `T` Air temperature (K) 
* `P` Air pressure (Pa) 
* `surface_type` Integer code for surface type (1=urban, 2=bare soil, etc.) 
* `vdep` Output: aerosol deposition velocity for surface (m/s) 




        

<hr>



### function canopy\_aero\_ddep\_subveg 

_Sub-canopy aerosol dry deposition velocity following Pleim et al._ 
```Fortran
subroutine, public canopy_aero_ddep_mod::canopy_aero_ddep_subveg (
    nlev nlev,
    z z,
    hc hc,
    lad lad,
    u u,
    d_p d_p,
    rho_p rho_p,
    T T,
    P P,
    vdep_opt vdep_opt,
    Ra Ra,
    modres modres,
    ustar ustar,
    vtype vtype,
    vdep vdep
) 
```



(2022),


Katul et al. (2010), Zhang et al. (2001), and Petroff et al. (2008) 

**Parameters:**


* `nlev` Number of canopy layers 
* `lad` Leaf area density profile (m^2/m^3) 
* `z` canopy model level heights (m) 
* `hc` canopy top level heights (m) 
* `u` Array of wind speed profile (m/s) 
* `d_p` Aerosol particle diameter (m) 
* `rho_p` Particle density (kg/m^3) 
* `T` Air temperature profile (K) 
* `P` Air pressure profile (Pa) 
* `vdep` Output: aerosol deposition velocity profile (m/s) 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_aero_ddep_mod.F90`

