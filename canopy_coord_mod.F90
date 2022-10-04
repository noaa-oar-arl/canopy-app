MODULE canopy_coord_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Coordinate and Domain Descriptions
! Purpose:  Contains canopy coordinate and domain descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------

    IMPLICIT NONE

!! .... defines geographic domain of inputs (read from file and/or user namelist)
    INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
    integer            :: nlat        !length of x coordinate
    integer            :: nlon        !length of y coordinate
    integer            :: modlays     !Number of total above and below canopy model layers
    real(rk)           :: modres      !Real value of model above and below canopy vertical resolution (m)

END MODULE canopy_coord_mod
