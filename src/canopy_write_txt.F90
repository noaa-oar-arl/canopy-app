
!> \file canopy_write_txt.F90
!! \brief Text Output File Writer Subroutine
!! \details This file contains the subroutine for writing canopy model output
!! to text format files. It serves as a wrapper that calls the main text
!! writing routine with appropriate parameters.
!!
!! \author Patrick C. Campbell
!! \date October 2022

!> \defgroup write_txt Text Output Writing
!! \brief Routines for writing text format output files
!! \{

!> \brief Write canopy outputs to text files
!! \details This subroutine writes canopy model calculation results to text format
!! output files. It serves as a wrapper for the main text writing routine,
!! passing the output file prefix and current time information.
!!
!! \param[in] OUTPREFX Output file prefix string
!! \param[in] TIMENOW Current time stamp string
SUBROUTINE canopy_write_txt(OUTPREFX,TIMENOW)

    USE canopy_txt_io_mod  !> main IO text reader/writer

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )  :: OUTPREFX  !> Output file prefix string
    CHARACTER(LEN=*), INTENT( IN )  :: TIMENOW   !> Current time stamp string

    !> \brief Call main text writing routine
    call write_txt(OUTPREFX,TIMENOW)

!> \}

END SUBROUTINE canopy_write_txt
