SUBROUTINE canopy_close_files(OUTPREFX)

!-------------------------------------------------------------------------------
! Name:     Close Files
! Purpose:  Close I/O API files.
! Revised:  10 Sep 2001  Original version.  (T. Otte)
!           09 Jan 2002  Changed "stop" statements to calls to "m3exit" for
!                        graceful shut-down of I/O API files.  (T. Otte)
!           29 Aug 2011  Replaced module IODECL3 with I/O API module M3UTILIO.
!                        Replaced F77 character declarations with F90 standard.
!                        Improved error handling.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           19 Dec 2018  Added runtime option to choose output format.
!                        (T. Spero)
!           04 Aug 2023  Using for canopy-app closing output files (P. Campbell)
!-------------------------------------------------------------------------------

    use canopy_files_mod
    use netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT( IN )    :: OUTPREFX
    CHARACTER(LEN=256)                :: fl
    CHARACTER(LEN=16),  PARAMETER     :: pname      = 'CLOSE_FILES'
    INTEGER                           :: rcode

!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

    CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CLOSING NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, '***   ', a, &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Gracefully close output files.
!-------------------------------------------------------------------------------

    fl = TRIM(OUTPREFX)//trim('.nc')

    rcode = nf90_close (cdfid_m)
    IF ( rcode /= nf90_noerr ) THEN
        WRITE (6,f9100) TRIM(pname), TRIM(fl),  &
            TRIM(nf90_strerror(rcode))
        CALL exit (2)
    ENDIF

END SUBROUTINE canopy_close_files
