

# File canopy\_write\_txt.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_write\_txt.F90**](canopy__write__txt_8F90.md)

[Go to the documentation of this file](canopy__write__txt_8F90.md)


```Fortran



SUBROUTINE canopy_write_txt(OUTPREFX,TIMENOW)

    USE canopy_txt_io_mod  !> main io text reader/writer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )  :: OUTPREFX
    CHARACTER(LEN=*), INTENT( IN )  :: TIMENOW

    call write_txt(outprefx,timenow)


END SUBROUTINE canopy_write_txt
```


