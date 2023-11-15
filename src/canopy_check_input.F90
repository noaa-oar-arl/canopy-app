
SUBROUTINE canopy_check_input(INFILE)

!-------------------------------------------------------------------------------
! Name:     Check Canopy Inputs (TXT or NETCDF)
! Purpose:  Check Canopy Inputs (TXT or NETCDF)
! Revised:  21 Dec 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_canopts_mod !main canopy option descriptions
    use canopy_ncf_io_mod, only: canopy_read_ncf

    implicit none


    !Local variables
    integer ppos
    CHARACTER(LEN=*), INTENT( IN )  :: INFILE

    ppos = scan(trim(INFILE),".", BACK= .true.)
    if (trim(INFILE(ppos:)).eq.".txt") then !TXT File
        if(infmt_opt .ne. 1) then !check to make sure input format matches text
            write(*,*)  'Wrong choice of INFMT_OPT ', infmt_opt, ' in namelist...exiting'
            write(*,*)  'Reading .txt file, change to INFMT_OPT = 1 '
            call exit(2)
        else !read text file
            call canopy_read_txt(INFILE)
            !TODO: if infmt_opt = 1 and i and j = 1 and var3d_opt=1 and pavd_opt=1 then...
            !call canopy_read_txt_3d(INFILE)...supporting 3D text file in canopy (e.g., PAVD)
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
