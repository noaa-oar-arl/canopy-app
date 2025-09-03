

# File canopy\_files\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_files\_mod.F90**](canopy__files__mod_8F90.md)

[Go to the documentation of this file](canopy__files__mod_8F90.md)


```Fortran



MODULE canopy_files_mod

    IMPLICIT NONE


    INTEGER                       :: cdfid_m
    INTEGER,            PARAMETER :: max_mm       = 10000
    INTEGER,            PARAMETER :: iutnml       =  8



    CHARACTER(LEN=256)            :: file_vars    ( max_mm )
    CHARACTER(LEN=256)            :: file_canvars ( max_mm )
    CHARACTER(LEN=256)            :: file_out     ( 1 )



    CHARACTER(LEN=*), PARAMETER   :: file_nml     = 'input/namelist.canopy'



END MODULE canopy_files_mod
```


