MODULE canopy_txt_io_mod

!-------------------------------------------------------------------------------
! Name:     TXT IO
! Purpose:  Contains routines to read met/sfc model text output.
! Revised:  03 Oct 2022  Original version  (P.C. Campbell)
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE read_met_txt(TXTFILE)

        USE canopy_coord_mod
        USE canopy_canmet_mod

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT( IN )  :: TXTFILE

        !Local variables
        integer i0, loc

        ! ... read met/sfc input variables from text file
        open(8,  file=TXTFILE,  status='old')
        i0 = 0
        read(8,*,iostat=i0)  ! skip headline
        do loc=1, nlat*nlon
            read(8, *) variables(loc)
        end do
        close(8)

    END SUBROUTINE read_met_txt

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

END MODULE canopy_txt_io_mod
