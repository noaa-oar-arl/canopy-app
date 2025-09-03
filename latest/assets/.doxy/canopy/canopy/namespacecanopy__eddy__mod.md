

# Namespace canopy\_eddy\_mod



[**Namespace List**](namespaces.md) **>** [**canopy\_eddy\_mod**](namespacecanopy__eddy__mod.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_eddyx**](#function-canopy_eddyx) (HCM HCM, ZK ZK, USTAR USTAR, MOL MOL, KZ KZ) <br>_Calculate canopy eddy diffusivity._  |




























## Public Functions Documentation




### function canopy\_eddyx 

_Calculate canopy eddy diffusivity._ 
```Fortran
subroutine canopy_eddy_mod::canopy_eddyx (
    HCM HCM,
    ZK ZK,
    USTAR USTAR,
    MOL MOL,
    KZ KZ
) 
```



Computes Eddy Diffusivity, Kz, within and just above canopy using stability-dependent parameterizations for different atmospheric conditions




**Parameters:**


* `HCM` Height of canopy top (m) 
* `ZK` Above/Below canopy height, z (m) 
* `USTAR` Model input friction velocity (m/s) 
* `MOL` Model input Monin-Obukhov Length (m) 
* `KZ` Estimated Eddy Diffusivity with canopy turbulence (m²/s) 



**Author:**

P. C. Campbell 




**Date:**

Jun 2022 




**Note:**

Based on Makar et al. (2017) algorithms for canopy turbulence Makar, P., Staebler, R., Akingunola, A. et al. (2017). The effects of forest canopy shading and turbulence on boundary layer ozone. Nature Communications, 8, 15243. [https://doi.org/10.1038/ncomms15243](https://doi.org/10.1038/ncomms15243) Raupach, M. R. (1989). A Practical Lagrangian method for relating scalar concentrations to source distributions in vegetation canopies. Quarterly Journal of the Royal Meteorological Society, 115, 609-632. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `src/canopy_eddy_mod.F90`

