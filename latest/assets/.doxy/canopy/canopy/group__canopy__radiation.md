

# Group canopy\_radiation



[**Modules**](modules.md) **>** [**canopy\_radiation**](group__canopy__radiation.md)



_Radiation attenuation and PPFD calculations within canopy._ 






































## Public Functions

| Type | Name |
| ---: | :--- |
|  subroutine | [**canopy\_rad\_mod::canopy\_fsun\_clu**](#function-canopy_fsun_clu) (FCLAI FCLAI, LAI LAI, CLU CLU, COSZEN COSZEN, FSUN FSUN) <br>_Calculate sunlit fraction using clumping index._  |
|  subroutine | [**canopy\_rad\_mod::canopy\_ppfd\_exp**](#function-canopy_ppfd_exp) (ZK ZK, FCH FCH, SFCRAD SFCRAD, LAI LAI, FSUN FSUN, PPFD\_SUN PPFD\_SUN, PPFD\_SHADE PPFD\_SHADE, PPFD\_AVE PPFD\_AVE) <br>_Calculate PPFD profiles using exponential model._  |




























## Public Functions Documentation




### function canopy\_fsun\_clu 

_Calculate sunlit fraction using clumping index._ 
```Fortran
subroutine canopy_rad_mod::canopy_fsun_clu (
    FCLAI FCLAI,
    LAI LAI,
    CLU CLU,
    COSZEN COSZEN,
    FSUN FSUN
) 
```



Computes sunlit/shaded fraction through canopy using photolysis correction factor and clumping index




**Parameters:**


* `FCLAI` Fractional cumulative LAI profile (dimensionless) 
* `LAI` Model input total Leaf Area Index (m²/m²) 
* `CLU` Model input Clumping Index (dimensionless) 
* `COSZEN` Model input Cosine Solar Zenith Angle (dimensionless) 
* `FSUN` Sunlit/Shaded fraction from photolysis correction factor 



**Author:**

P. C. Campbell 




**Date:**

Jun 2023 




**Note:**

Based on Bonan (2019) equation 14.18 for clumping correction 





        

<hr>



### function canopy\_ppfd\_exp 

_Calculate PPFD profiles using exponential model._ 
```Fortran
subroutine canopy_rad_mod::canopy_ppfd_exp (
    ZK ZK,
    FCH FCH,
    SFCRAD SFCRAD,
    LAI LAI,
    FSUN FSUN,
    PPFD_SUN PPFD_SUN,
    PPFD_SHADE PPFD_SHADE,
    PPFD_AVE PPFD_AVE
) 
```



Computes photosynthetic photon flux density for sunlit and shaded leaves through canopy using exponential decay models from Silva et al. (2020)




**Parameters:**


* `ZK` Input model heights (m) 
* `FCH` Model input canopy height (m) 
* `SFCRAD` Model input instantaneous surface downward shortwave flux (W/m²) 
* `LAI` Model input total Leaf Area Index (m²/m²) 
* `FSUN` Sunlit/Shaded fraction from photolysis correction factor 
* `PPFD_SUN` PPFD for sunlit leaves (μmol photons/m²/s) 
* `PPFD_SHADE` PPFD for shaded leaves (μmol photons/m²/s) 
* `PPFD_AVE` Average PPFD for sunlit and shaded leaves (μmol photons/m²/s) 



**Author:**

P. C. Campbell 




**Date:**

Jun 2023 




**Note:**

Based on Silva et al. (2020) 5-layer canopy exponential PPFD model Silva, S. J., Heald, C. L., and Guenther, A. B. (2020). Development of a reduced-complexity plant canopy physics surrogate model for use in chemical transport models: a case study with GEOS-Chem v12.3.0. Geoscientific Model Development, 13, 2569-2585. [https://doi.org/10.5194/gmd-13-2569-2020](https://doi.org/10.5194/gmd-13-2569-2020) 





        

<hr>

------------------------------


