
SUBROUTINE canopy_write_txt(OUTPREFX,TIMENOW)

!-------------------------------------------------------------------------------
! Name:     Write Canopy Outputs to TXT
! Purpose:  Write Canopy Outputs to TXT
! Revised:  06 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_txt_io_mod  !main IO text reader/writer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )  :: OUTPREFX
    CHARACTER(LEN=*), INTENT( IN )  :: TIMENOW

    call write_txt(OUTPREFX,TIMENOW)


END SUBROUTINE canopy_write_txt
