
SUBROUTINE canopy_readnml (nlat,nlon,canlays,canres,href,z0ghcm,lamdars)

!-------------------------------------------------------------------------------
! Name:     Read Canopy Namelist
! Purpose:  Reads input namelist to get user control variables.
!           15 Jul 2022  Original Version (P.C. Campbell)
!                        
!-------------------------------------------------------------------------------

  USE canopy_files_mod

  IMPLICIT NONE

  INTEGER,               INTENT(OUT) :: nlat,nlon,canlays
  REAL,                  INTENT(OUT) :: canres,href,z0ghcm,lamdars
  INTEGER                            :: istat
  INTEGER                            :: n
  CHARACTER(LEN=16),     PARAMETER   :: pname      = 'CANOPY_READNML'

  NAMELIST /filenames/   file_prof, file_vars

  NAMELIST /userdefs/    nlat, nlon, canlays, canres, href, z0ghcm, lamdars

!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING CANOPY NAMELIST FILE ON UNIT ', i3, &
    & /, 1x, '***   NAMELIST FILE NAME = ', a, &
    & /, 1x, '***   IOSTAT = ', i4, &
    & /, 1x, 70('*'))"


  CHARACTER(LEN=256), PARAMETER :: f9050 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR READING NAMELIST FILE ON UNIT ', i3, &
    & /, 1x, '***   NAMELIST FILE NAME = ', a, &
    & /, 1x, '***   NAMELIST = ', a, &
    & /, 1x, '***   IOSTAT = ', i4, &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Open canopy namelist file.
!-------------------------------------------------------------------------------

  OPEN (iutnml, FILE=file_nml, STATUS='OLD', IOSTAT=istat)

  IF ( istat > 0 ) THEN
    WRITE (*,f9000) TRIM(pname), iutnml, TRIM(file_nml), istat
    CALL EXIT(istat)
  ENDIF

!-------------------------------------------------------------------------------
! Initialize canopy input file names.
!-------------------------------------------------------------------------------

  file_prof(:)    = " "
  file_vars(:)    = " "

!-------------------------------------------------------------------------------

! Set default value for number of latitude/longitude cells (default = 1D point)
  nlat = 1
  nlon = 1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer value for number of canopy layers (default = 100 layers)
  canlays = 100
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for canopy vertical resolution (m) (Default = 0.5 m)
  canres = 0.5
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for reference height above canopy (m) (Default = 10 m)
  href  = 10.0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for ratio of ground roughness length to canopy top height
  z0ghcm  = 0.0025
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for Influence function associated with roughness sublayer
  lamdars  = 1.25
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read namelist to get user definitions.  Rewind namelist file after each
! read in case namelists are not in the correct order in the namelist.
!-------------------------------------------------------------------------------

  READ (iutnml, filenames, IOSTAT=istat)
  IF ( istat > 0 ) THEN
    WRITE (*,f9050) TRIM(pname), iutnml, TRIM(file_nml), "filenames", istat
    CALL EXIT(istat)
  ENDIF
  REWIND (iutnml)

  READ (iutnml, userdefs, IOSTAT=istat)
  IF ( istat > 0 ) THEN
    WRITE (*,f9050) TRIM(pname), iutnml, TRIM(file_nml), "userdefs", istat
    CALL EXIT(istat)
  ENDIF
  REWIND (iutnml)

!-------------------------------------------------------------------------------
! Crop blank spaces off ends of file names.
!-------------------------------------------------------------------------------

  DO n = 1, SIZE(file_prof)
    file_prof(n) = TRIM( ADJUSTL( file_prof(n) ) )
    file_vars(n)= TRIM( ADJUSTL( file_vars(n) ) )
  ENDDO

!-------------------------------------------------------------------------------
! Verify values of user-defined options (need conditions added...)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Close namelist file.
!-------------------------------------------------------------------------------
!  write(*,*)'namelist intvl=',intvl
  CLOSE (iutnml)

END SUBROUTINE canopy_readnml
