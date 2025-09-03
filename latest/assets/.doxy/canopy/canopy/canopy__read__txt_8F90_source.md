

# File canopy\_read\_txt.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_read\_txt.F90**](canopy__read__txt_8F90.md)

[Go to the documentation of this file](canopy__read__txt_8F90.md)


```Fortran



SUBROUTINE canopy_read_txt(INFILE,INFILE2)

    USE canopy_canopts_mod !> main canopy options
    USE canopy_txt_io_mod  !> main io text reader/writer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )  :: INFILE
    CHARACTER(LEN=*), INTENT( IN )  :: INFILE2

    if (var3d_opt .ne. 1) then ! reading only variable text file
        call read_txt(infile)
    else                       ! reading variable text file and canopy profile
        call read_txt(infile)
        call read_can_txt(infile2)
    end if


END SUBROUTINE canopy_read_txt
```


