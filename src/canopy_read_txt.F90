
!> \file canopy_read_txt.F90
!! \brief Text Input File Reader Subroutine
!! \details This file contains the subroutine for reading meteorological and surface
!! data from text format input files. It handles both main variable files and
!! supplementary canopy profile files.
!!
!! \author Patrick C. Campbell
!! \date October 2022

!> \defgroup read_txt Text Input Reading
!! \brief Routines for reading text format input files
!! \{

!> \brief Read canopy met/sfc inputs from text files
!! \details This subroutine reads meteorological and surface input data from text files.
!! It can handle both single variable text files and combined variable plus canopy
!! profile text files depending on the var3d_opt setting.
!!
!! \param[in] INFILE Primary input text file path
!! \param[in] INFILE2 Secondary canopy profile input text file path
SUBROUTINE canopy_read_txt(INFILE,INFILE2)

    USE canopy_canopts_mod !> main canopy options
    USE canopy_txt_io_mod  !> main IO text reader/writer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )  :: INFILE   !> Primary input text file path
    CHARACTER(LEN=*), INTENT( IN )  :: INFILE2  !> Secondary canopy profile input text file path

    !> \brief Determine which text files to read based on 3D variable option
    if (var3d_opt .ne. 1) then ! reading only variable text file
        call read_txt(INFILE)
    else                       ! reading variable text file and canopy profile
        call read_txt(INFILE)
        call read_can_txt(INFILE2)
    end if

!> \}

END SUBROUTINE canopy_read_txt
