
SUBROUTINE canopy_read_txt(INFILE)

!-------------------------------------------------------------------------------
! Name:     Read Canopy Met/Sfc Inputs from TXT
! Purpose:  Read Canopy Met/Sfc Inputs from TXT
! Revised:  03 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_txt_io_mod  !main IO text reader/writer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )  :: INFILE


    call read_txt(INFILE)


END SUBROUTINE canopy_read_txt
