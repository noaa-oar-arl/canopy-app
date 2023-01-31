SUBROUTINE canopy_readnml

!-------------------------------------------------------------------------------
! Name:     Read Canopy Namelist
! Purpose:  Reads input namelist to get user control variables.
!           15 Jul 2022  Original Version (P.C. Campbell)
!
!-------------------------------------------------------------------------------
    USE canopy_files_mod
    USE canopy_canopts_mod
    USE canopy_coord_mod
    use canopy_canvars_mod, ONLY: zk

    IMPLICIT NONE

    !local variables
    INTEGER                            :: istat
    INTEGER                            :: n,i
    CHARACTER(LEN=*),      PARAMETER   :: pname = 'CANOPY_READNML'

    NAMELIST /filenames/ file_vars, file_out

    NAMELIST /userdefs/  infmt_opt, nlat, nlon, modlays, modres, href_opt, href_set, &
        z0ghc, lambdars, flameh_opt, flameh_set, ifcanwind, ifcanwaf, ifcaneddy, &
        ifcanphot, pai_opt, pai_set, lu_opt, z0_opt, dx_opt, dx_set, lai_thresh, &
        frt_thresh, fch_thresh, rsl_opt

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
! Initialize canopy file names.
!-------------------------------------------------------------------------------

    file_vars(:) = " "
    file_out(:)  = " "

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for 1D (txt or ncf) or 2D (ncf only) input file format
! (default = 0, i.e., 2D)
    infmt_opt = 0
!-------------------------------------------------------------------------------

! Set default value for number of latitude/longitude cells (default = 1D point)
    nlat = 1
    nlon = 1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer value for number of canopy layers (default = 100 layers)
    modlays = 100
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for canopy vertical resolution (m) (Default = 0.5 m)
    modres = 0.5_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for reference height set values or array from file (default = 0)
    href_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default value for reference height above canopy (m) (Default = 10 m)
    href_set = 10.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for ratio of ground roughness length to canopy top height
    z0ghc = 0.0025_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for Influence function associated with roughness sublayer
    lambdars = 1.25_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for flame height set values or calculation (default = 0)
    flameh_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for flame height (m) (Default = 2.0 m)
    flameh_set = 2.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy wind option (default = .FALSE.)
    ifcanwind = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy WAF option (default = .FALSE.)
    ifcanwaf = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy eddy diffusivity (default = .FALSE.)
    ifcaneddy = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default logical for canopy photolysis attenuation (default = .FALSE.)
    ifcanphot = .FALSE.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for PAI set values or calculation (default = 0)
    pai_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for PAI set value (default = 4.0)
    pai_set = 4.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set integer for LU type used from model mapped to Massman et al. (default = 0)
    lu_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set integer for z0 estimate either from model or vegtype dependent (default = 0)
    z0_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default integer for DX set values or calculation (default = 0)
    dx_opt = 0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for DX Cell resolution set value (default = 1.0 m)
    dx_set = 1.0_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for user LAI threshold value for canopy (default = 0.1)
    lai_thresh = 0.1_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for user FRT threshold value for canopy (default = 0.5)
    frt_thresh = 0.5_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set default real value for user FCH threshold value for canopy (default = 0.5 m)
    fch_thresh = 0.1_rk
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Set integer for unified RSL option used in model (default = 0, off)
    rsl_opt = 0
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

    DO n = 1, SIZE(file_vars)
        file_vars(n)= TRIM( ADJUSTL( file_vars(n) ) )
    ENDDO

    DO n = 1, SIZE(file_out)
        file_out(n)= TRIM( ADJUSTL( file_out(n) ) )
    ENDDO

!-------------------------------------------------------------------------------
! Verify values of user-defined options (need conditions added...)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Derive canopy model profile heights from user NL input
!-------------------------------------------------------------------------------

    if(.not.allocated(zk))         allocate(zk(modlays))
    zk(1) = 0.0
    do i=2, modlays
        zk(i)   = zk(i-1) + modres
    end do

!-------------------------------------------------------------------------------
! Close namelist file.
!-------------------------------------------------------------------------------
!  write(*,*)'namelist intvl=',intvl
    CLOSE (iutnml)

END SUBROUTINE canopy_readnml
