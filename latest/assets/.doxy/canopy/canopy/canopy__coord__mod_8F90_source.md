

# File canopy\_coord\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_coord\_mod.F90**](canopy__coord__mod_8F90.md)

[Go to the documentation of this file](canopy__coord__mod_8F90.md)


```Fortran



MODULE canopy_coord_mod
    use canopy_const_mod, ONLY: rk
    IMPLICIT NONE


    CHARACTER(LEN=24)  :: time_start
    CHARACTER(LEN=24)  :: time_end
    INTEGER            :: time_intvl
    integer            :: ntime



    integer            :: nlat
    integer            :: nlon



    integer            :: modlays
    real(rk)           :: modres



END MODULE canopy_coord_mod
```


