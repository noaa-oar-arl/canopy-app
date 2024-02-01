
SUBROUTINE canopy_check_input(INFILE,INFILE2)

!-------------------------------------------------------------------------------
! Name:     Check Canopy Inputs (TXT or NETCDF)
! Purpose:  Check Canopy Inputs (TXT or NETCDF)
! Revised:  21 Dec 2022  Original version.  (P.C. Campbell)
! Revised:  30 Nov 2023  Added supplementary canopy profile, INFILE2.  (P.C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_canopts_mod !main canopy option descriptions
    use canopy_ncf_io_mod, only: canopy_read_ncf

    implicit none


    !Local variables
    integer ppos
    CHARACTER(LEN=*), INTENT( IN )  :: INFILE,INFILE2

    ppos = scan(trim(INFILE),".", BACK= .true.)
    if (trim(INFILE(ppos:)).eq.".txt") then !TXT File
        if(infmt_opt .ne. 1) then !check to make sure input format matches text
            write(*,*)  'Wrong choice of INFMT_OPT ', infmt_opt, ' in namelist...exiting'
            write(*,*)  'Reading .txt file, change to INFMT_OPT = 1 '
            call exit(2)
        else !read text file
            call canopy_read_txt(INFILE,INFILE2)
        end if
    else if (trim(INFILE(ppos:)).eq.".nc") then !NetCDF File
        call canopy_read_ncf(INFILE)
    else if (trim(INFILE(ppos:)).eq.".ncf") then
        call canopy_read_ncf(INFILE)
    else if (trim(INFILE(ppos:)).eq.".nc4") then
        call canopy_read_ncf(INFILE)
    else
        write(*,*)  'Error the file input type ',trim(INFILE(ppos:)), &
            ' is not supported...exiting'
        call exit(2)
    end if   !File Input types

END SUBROUTINE canopy_check_input
