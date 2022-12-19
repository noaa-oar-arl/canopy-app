
SUBROUTINE canopy_read_ncf(infile)

!-------------------------------------------------------------------------------
! Name:     Read Canopy Met/Sfc Inputs from NetCDF
! Purpose:  Read Canopy Met/Sfc Inputs from NetCDF
! Revised:  19 Dec 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_ncf_io_mod  !main IO NetCDF reader
    USE netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )   :: infile
    CHARACTER(LEN=*),  PARAMETER     :: pname      = 'canopy_read_ncf'
    INTEGER                          :: cdfid, rcode

!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

    CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM NCF INPUT FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   RCODE = ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9410 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING NCF ID FROM NCF INPUT FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9420 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR INQUIRING ABOUT VAR IN NCF INPUT FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9430 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING DIMS FROM NCF INPUT FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

    CHARACTER(LEN=256), PARAMETER :: f9440 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING NCF INPUT FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"


!-------------------------------------------------------------------------------
! Open NCF input file
!-------------------------------------------------------------------------------

    rcode = nf90_open (infile, nf90_nowrite, cdfid)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9440) TRIM(pname), TRIM(infile)
        call exit(2)
    ENDIF

!-------------------------------------------------------------------------------
! Extract and check input lat/lon dimensions and time information.
!-------------------------------------------------------------------------------


END SUBROUTINE canopy_read_ncf
