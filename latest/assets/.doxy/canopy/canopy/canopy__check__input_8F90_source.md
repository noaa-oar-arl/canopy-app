

# File canopy\_check\_input.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_check\_input.F90**](canopy__check__input_8F90.md)

[Go to the documentation of this file](canopy__check__input_8F90.md)


```Fortran



SUBROUTINE canopy_check_input(INFILE,INFILE2)
    use canopy_canopts_mod !> main canopy option descriptions
    use canopy_ncf_io_mod, only: canopy_read_ncf

    implicit none

    integer ppos
    CHARACTER(LEN=*), INTENT( IN )  :: infile
    CHARACTER(LEN=*), INTENT( IN )  :: infile2

    ppos = scan(trim(infile),".", back= .true.)
    if (trim(infile(ppos:)).eq.".txt") then !TXT File
        if(infmt_opt .ne. 1) then !check to make sure input format matches text
            write(*,*)  'Wrong choice of INFMT_OPT ', infmt_opt, ' in namelist...exiting'
            write(*,*)  'Reading .txt file, change to INFMT_OPT = 1 '
            call exit(2)
        else !read text file
            call canopy_read_txt(infile,infile2)
        end if
    else if (trim(infile(ppos:)).eq.".nc") then !NetCDF File
        call canopy_read_ncf(infile)
    else if (trim(infile(ppos:)).eq.".ncf") then
        call canopy_read_ncf(infile)
    else if (trim(infile(ppos:)).eq.".nc4") then
        call canopy_read_ncf(infile)
    else
        write(*,*)  'Error the file input type ',trim(infile(ppos:)), &
            ' is not supported...exiting'
        call exit(2)
    end if   !File Input types


END SUBROUTINE canopy_check_input
```


