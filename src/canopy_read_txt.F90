
SUBROUTINE canopy_read_txt(INFILE,INFILE2)

!-------------------------------------------------------------------------------
! Name:     Read Canopy Met/Sfc Inputs from TXT
! Purpose:  Read Canopy Met/Sfc Inputs from TXT
! Revised:  03 Oct 2022  Original version.  (P.C. Campbell)
! Revised:  30 Nov 2023  Added supplementary canopy profile.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_canopts_mod !main canopy options
    USE canopy_txt_io_mod  !main IO text reader/writer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )  :: INFILE,INFILE2

    if (var3d_opt .ne. 1) then ! reading only variable text file
        call read_txt(INFILE)
    else                       ! reading variable text file and canopy profile
        call read_txt(INFILE)
        call read_can_txt(INFILE2)
    end if

END SUBROUTINE canopy_read_txt
