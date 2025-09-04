
!> \file canopy_check_input.F90
!! \brief Input File Validation and Reading Subroutine
!! \details This file contains the subroutine that checks and reads canopy input files.
!! It supports both text and NetCDF input formats and validates the file format
!! against the user-specified input format option.
!!
!! \author Patrick C. Campbell
!! \date December 2022

!> \defgroup check_input Input File Validation
!! \brief Routines for validating and reading input files
!! \{

!> \brief Check and read canopy input files (TXT or NETCDF)
!! \details This subroutine determines the input file format based on the file extension
!! and validates it against the user-specified format option (infmt_opt). It then
!! calls the appropriate reading routine:
!! - For .txt files: calls canopy_read_txt()
!! - For .nc/.ncf/.nc4 files: calls canopy_read_ncf()
!!
!! The subroutine performs error checking to ensure the file format matches the
!! namelist specification and exits with an error code if there are mismatches.
!!
!! \param[in] INFILE Primary input file path
!! \param[in] INFILE2 Secondary canopy profile input file path
SUBROUTINE canopy_check_input(INFILE,INFILE2)
    use canopy_canopts_mod !> main canopy option descriptions
    use canopy_ncf_io_mod, only: canopy_read_ncf

    implicit none

!> \defgroup check_input_vars Local Variables
!! \brief Local variables for input file processing
!! \{
    integer ppos                                !> Position of file extension separator
    CHARACTER(LEN=*), INTENT( IN )  :: INFILE  !> Primary input file path
    CHARACTER(LEN=*), INTENT( IN )  :: INFILE2 !> Secondary canopy profile input file path
!> \}

    !> \brief Determine file format and validate against namelist option
    !! \details Find the file extension and check if it matches the user-specified format option
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

!> \}

END SUBROUTINE canopy_check_input
