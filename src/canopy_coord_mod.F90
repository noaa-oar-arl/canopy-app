MODULE canopy_coord_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Coordinate and Domain Descriptions
! Purpose:  Contains canopy coordinate and domain descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!           03 Aug 2023  Added ntim (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk
    IMPLICIT NONE

!! .... defines time and geographic domain of inputs (read from file and/or user namelist)
    CHARACTER(LEN=24)  :: time_start  !user input: YYYY-MM-DD-HH:MM:SS.SSSS
    CHARACTER(LEN=24)  :: time_end    !user input: YYYY-MM-DD-HH:MM:SS.SSSS
    INTEGER            :: time_intvl  !user input: time interval for input/output [hrs]
    integer            :: ntime       !number of model timesteps
    integer            :: nlat        !length of x coordinate
    integer            :: nlon        !length of y coordinate
    integer            :: modlays     !Number of total above and below canopy model layers
    real(rk)           :: modres      !Real value of model above and below canopy vertical resolution (m)

END MODULE canopy_coord_mod
