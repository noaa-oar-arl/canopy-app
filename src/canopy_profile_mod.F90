module canopy_profile_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_PARM( VTYPE, FCH, FFRAC, LAI, &
        PAI_OPT, PAI_SET, LU_OPT, FIRETYPE, CDRAG, &
        PAI, ZCANMAX, SIGMAU, SIGMA1 )

!-----------------------------------------------------------------------

! Description:
!     determines the canopy, fire, and foliage shape distribution parameters

! Preconditions:
!     lat, lon, and vegetation type

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial canopy parameter model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

        use canopy_const_mod, ONLY: rk      !constants for canopy models
        use canopy_utils_mod, ONLY: CalcPAI !canopy utilities/functions

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )  :: VTYPE           ! Grid cell dominant vegetation type
        REAL(RK),    INTENT( IN )  :: FCH             ! Grid cell canopy height (m)
        REAL(RK),    INTENT( IN )  :: FFRAC           ! Grid cell forest fraction
        REAL(RK),    INTENT( IN )  :: LAI             ! Grid cell leaf area index
        INTEGER,     INTENT( IN )  :: PAI_OPT         ! integer for PAI values used or calculated (default = 0)
        REAL(RK),    INTENT( IN )  :: PAI_SET         ! real value for PAI set values used (default = 4.0)
        INTEGER,     INTENT( IN )  :: LU_OPT          ! integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)

        INTEGER,     INTENT( OUT ) :: FIRETYPE        ! 1 = Above Canopy Fire; 0 = Below Canopy Fire; -1 No Canopy
        REAL(RK),    INTENT( OUT ) :: CDRAG           ! Drag coefficient (nondimensional)
        REAL(RK),    INTENT( OUT ) :: PAI             ! Plant/foliage area index (nondimensional)
        REAL(RK),    INTENT( OUT ) :: ZCANMAX         ! Height of maximum foliage area density (z/h) (nondimensional)
        REAL(RK),    INTENT( OUT ) :: SIGMAU          ! Standard deviation of shape function above zcanmax (z/h)
        REAL(RK),    INTENT( OUT ) :: SIGMA1          ! Standard deviation of shape function below zcanmax (z/h)

        if (LU_OPT .eq. 0 .or. LU_OPT .eq. 1) then !VIIRS or MODIS LU types

            !approx/average vegtype mapping to Massman et al. forest types
            if (VTYPE .ge. 1 .and. VTYPE .le. 2) then !VIIRS/MODIS Cat 1-2/Evergreen Needleleaf & Broadleaf
                !--> Use average Massman Spruce+ScotsPine+JackPine+LoblollyPine Forest
                FIRETYPE=0
                CDRAG=(0.25_rk + 0.20_rk + 0.20_rk + 0.20_rk)/4.0_rk
                if (PAI_OPT .eq. 0) then      !Katul et al. 2004 vegtype
                    PAI=(3.28_rk + 2.41_rk + 2.14_rk + 3.78_rk)/4.0_rk
                else if (PAI_OPT .eq. 1) then !PAI calculation (Massman et al., Eq. 19)
                    PAI=CalcPAI(FCH,FFRAC)
                else if (PAI_OPT .eq. 2) then !PAI=LAI/(1-alpha), where alpha is the "woody-to-total area ratio"
                    !and is vegetation type dependent from Fang et al. (2019),
                    !https://doi.org/10.1029/2018RG000608:
                    PAI=(LAI/(1.0_rk - 0.2_rk)) !Assume alpha evergreen are more boreal/conifer softwoods
                else if (PAI_OPT .eq. 3) then !PAI value from user
                    PAI=PAI_SET
                else
                    write(*,*)  'Wrong PAI_OPT choice of ', PAI_OPT, 'in namelist...exiting'
                    call exit(2)
                end if
                ZCANMAX=(0.60_rk + 0.60_rk + 0.58_rk + 0.60_rk)/4.0_rk
                SIGMAU=(0.38_rk  + 0.30_rk + 0.20_rk + 0.10_rk)/4.0_rk
                SIGMA1=(0.16_rk  + 0.10_rk + 0.20_rk + 0.27_rk)/4.0_rk
            end if

            if (VTYPE .ge. 3 .and. VTYPE .le. 4) then !VIIRS/MODIS Cat 3-4 Deciduous Needleleaf and  Broadleaf
                !--> Use Massman Hardwood Forest + Aspen
                FIRETYPE=0
                CDRAG=(0.15_rk + 0.20_rk)/2.0_rk
                if (PAI_OPT .eq. 0) then      !Katul et al. 2004 vegtype
                    PAI=(4.93_rk + 3.28)/2.0_rk
                else if (PAI_OPT .eq. 1) then !Massman PAI calculation (Eq. 19)
                    PAI=CalcPAI(FCH,FFRAC)
                else if (PAI_OPT .eq. 2) then !PAI=LAI/(1-alpha), where alpha is the "woody-to-total area ratio"
                    !and is vegetation type dependent from Fang et al. (2019),
                    !https://doi.org/10.1029/2018RG000608:
                    PAI=(LAI/(1.0_rk - 0.1_rk)) !assume alpha deciduous have more tropical hardwoods
                else if (PAI_OPT .eq. 3) then !PAI value from user
                    PAI=PAI_SET
                else
                    write(*,*)  'Wrong PAI_OPT choice of ', PAI_OPT, 'in namelist...exiting'
                    call exit(2)
                end if
                ZCANMAX=(0.84_rk + 0.36_rk)/2.0_rk
                SIGMAU=(0.13_rk + 0.60_rk)/2.0_rk
                SIGMA1=(0.30_rk + 0.20_rk)/2.0_rk
            end if

            if (VTYPE .eq. 5) then !VIIRS/MODIS Cat 5 Mixed Forests
                !--> Use average Massman Aspen+Spruce+ScotsPine+JackPine+LoblollyPine+Hardwood Forest
                FIRETYPE=0
                CDRAG=(0.20_rk + 0.25_rk + 0.20_rk + 0.20_rk + 0.20_rk + 0.15_rk)/6.0_rk
                if (PAI_OPT .eq. 0) then      !Katul et al. 2004 vegtype
                    PAI=(5.73_rk + 3.28_rk + 2.41_rk + 2.14_rk + 3.78_rk + 4.93_rk)/6.0_rk
                else if (PAI_OPT .eq. 1) then !PAI calculation (Massman et al., Eq. 19)
                    PAI=CalcPAI(FCH,FFRAC)
                else if (PAI_OPT .eq. 2) then !PAI=LAI/(1-alpha), where alpha is the "woody-to-total area ratio"
                    !and is vegetation type dependent from Fang et al. (2019),
                    !https://doi.org/10.1029/2018RG000608:
                    PAI=(LAI/(1.0_rk - 0.15_rk))!assume alpha is avg. of evergreen and deciduous
                else if (PAI_OPT .eq. 3) then !PAI value from user
                    PAI=PAI_SET
                else
                    write(*,*)  'Wrong PAI_OPT choice of ', PAI_OPT, 'in namelist...exiting'
                    call exit(2)
                end if
                ZCANMAX=(0.60_rk + 0.36_rk + 0.60_rk + 0.58_rk + 0.60_rk + 0.84_rk)/6.0_rk
                SIGMAU=(0.38_rk + 0.60_rk + 0.30_rk + 0.20_rk + 0.10_rk + 0.13_rk)/6.0_rk
                SIGMA1=(0.16_rk + 0.20_rk + 0.10_rk + 0.20_rk + 0.27_rk + 0.30_rk)/6.0_rk
            end if

            if ((VTYPE .ge. 6 .and. VTYPE .le. 10) .or. VTYPE .eq. 12 ) then !VIIRS/MODIS Cat 6-10 or 12/Shrubs, Croplands, and Grasses
                !--> Average of Massman Corn + Rice )
                FIRETYPE=1
                CDRAG=(0.30_rk + 0.30_rk)/2.0_rk
                if (PAI_OPT .eq. 0) then      !Katul et al. 2004 vegtype
                    PAI=(2.94_rk + 3.10_rk)/2.0_rk
                else if (PAI_OPT .eq. 1) then !PAI calculation (Massman et al., Eq. 19)
                    PAI=CalcPAI(FCH,FFRAC)
                else if (PAI_OPT .eq. 2) then !PAI=LAI/(1-alpha), where alpha is the "woody-to-total area ratio"
                    !and is vegetation type dependent from Fang et al. (2019),
                    !https://doi.org/10.1029/2018RG000608:
                    PAI=(LAI/(1.0_rk - 0.32_rk))!assume alpha is avg. of shrubs + savanna
                else if (PAI_OPT .eq. 3) then !PAI value from user
                    PAI=PAI_SET
                else
                    write(*,*)  'Wrong PAI_OPT choice of ', PAI_OPT, 'in namelist...exiting'
                    call exit(2)
                end if
                ZCANMAX=(0.94_rk + 0.62_rk)/2.0_rk
                SIGMAU=(0.03_rk + 0.50_rk)/2.0_rk
                SIGMA1=(0.60_rk + 0.45_rk)/2.0_rk
            end if

        else
            write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS/MODIS available right now...exiting'
            call exit(2)
        end if


    END SUBROUTINE CANOPY_PARM
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_FOLIAGE( MODLAYS, ZHC, ZCANMAX, SIGMAU, SIGMA1, &
        FAFRACZINT )

!-----------------------------------------------------------------------

! Description:
!     computes canopy foliage and plant distribution functions.

! Preconditions:
!     height of canopy, and canopy distribution parameters

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Oct 2022 P.C. Campbell: Initial version
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk                   !constants for canopy models
        use canopy_utils_mod, ONLY: IntegrateTrapezoid   !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )  :: MODLAYS               ! Number of total above and below canopy model layers
        REAL(RK),    INTENT( IN )  :: ZCANMAX               ! Height of maximum foliage area density (z/h) (nondimensional)
        REAL(RK),    INTENT( IN )  :: SIGMAU                ! Standard deviation of shape function above zcanmax (z/h)
        REAL(RK),    INTENT( IN )  :: SIGMA1                ! Standard deviation of shape function below zcanmax (z/h)
        REAL(RK),    INTENT( IN )  :: ZHC(:)                ! z/h (dimensionless)
        REAL(RK),    INTENT( OUT ) :: FAFRACZINT(:)         ! integral of incremental fractional foliage shape function

!     Local variables
        INTEGER                     :: i
        REAL(RK), allocatable       :: fainc(:)              ! incremental foliage shape function
        REAL(RK), allocatable       :: fafracz(:)            ! incremental fractional foliage shape function
        REAL(RK)                    :: fatot                 ! integral of total fractional foliage shape function

! Citation:
!  W.J. Massman, J.M. Forthofer, and M.A. Finney. An improved
!  canopy wind model for predicting wind adjustment factors
!  and wildland fire behavior. Canadian Journal of Forest Research.
!  47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354

!Eqs. 1-2 on pg 595 of Massman et al., 2017

        if(.not.allocated(fainc))      allocate(fainc(MODLAYS))
        if(.not.allocated(fafracz))    allocate(fafracz(MODLAYS))

! ... initialize canopy profile dependent variables
        fainc             = 0.0_rk
        fafracz           = 0.0_rk

! ... calculate canopy/foliage distribution shape profile - bottom up total in-canopy and fraction at z
        do i=1, MODLAYS
            if (ZHC(i) >= ZCANMAX .and. ZHC(i) <= 1.0) then
                fainc(i) = exp((-1.0*((ZHC(i)-ZCANMAX)**2.0))/SIGMAU**2.0)
            else if (ZHC(i) >= 0.0 .and. ZHC(i) <= ZCANMAX) then
                fainc(i) = exp((-1.0*((ZCANMAX-ZHC(i))**2.0))/SIGMA1**2.0)
            end if
        end do

        fatot = IntegrateTrapezoid(ZHC,fainc)

! ... calculate plant distribution function
        do i=1, MODLAYS
            fafracz(i) = fainc(i)/fatot
            FAFRACZINT(i) = IntegrateTrapezoid(ZHC(1:i),fafracz(1:i))
        end do

    END SUBROUTINE CANOPY_FOLIAGE

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_ZPD( ZHC, FCLAI, UBZREF, Z0GHC, &
        LAMBDARS, CDRAG, PAI, FCH, HREF, Z0_MOD, &
        VTYPE, LU_OPT, Z0_OPT, d_h, zo_h )

!-----------------------------------------------------------------------

! Description:
!     computes zero plane displacement height and total surface roughness length basd on Massman et al. (2017).

! Preconditions:
!     foliage characteristics, drag, and in-canopy plant distribution functions

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone canopy wind model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk, vonk    !constants for canopy models
        use canopy_utils_mod     !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: ZHC     ( : )   ! SUB-CANOPY Total z/h layers (nondimensional)
        REAL(RK),    INTENT( IN )  :: FCLAI   ( : )   ! SUB-CANOPY Fractional (z) shapes of the
        ! plant surface distribution, i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )  :: UBZREF          ! Mean wind speed at reference height (ZREF) (m/s)
        REAL(RK),    INTENT( IN )  :: Z0GHC           ! Ratio of ground roughness length to canopy top height (nondimensional)
        REAL(RK),    INTENT( IN )  :: LAMBDARS         ! Value representing influence of roughness sublayer (nondimensional)
        REAL(RK),    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
        REAL(RK),    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
        REAL(RK),    INTENT( IN )  :: FCH             ! Grid cell canopy height (m)
        REAL(RK),    INTENT( IN )  :: HREF            ! Reference Height above the canopy (ZREF = FCH + HREF; m)
        REAL(RK),    INTENT( IN )  :: Z0_MOD          ! Input model value of surface roughness length, z0 (m)
        INTEGER,     INTENT( IN )  :: VTYPE           ! Grid cell dominant vegetation type
        INTEGER,     INTENT( IN )  :: LU_OPT          ! integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
        INTEGER,     INTENT( IN )  :: Z0_OPT          ! integer for setting first estimate of z0 (default = 0 for Z0_MOD)
        REAL(RK),    INTENT( OUT ) :: d_h             ! zero plane displacement height d/h (i.e., Product of dha*dhb)
        REAL(RK),    INTENT( OUT ) :: zo_h            ! Surface (soil+veg) roughness length, zo/h

!     Local variables
        real(rk)                   :: ustrmod         ! Friction Velocity parameterization (m/s)
        real(rk)                   :: dha             ! First term A in zero plane displacement height (d/h) (nondimensional)
        real(rk)                   :: dhb             ! Second term B in zero plane displacement height (d/h) (nondimensional)
        real(rk)                   :: qa              ! Coefficient term in Reynolds stress for d/h
        real(rk)                   :: qb              ! Coefficient term in Reynolds stress for d/h
        real(rk)                   :: qc              ! Coefficient term in Reynolds stress for d/h
        real(rk)                   :: qstar           ! Total Reynolds stress term
        real(rk)                   :: fafraczInt_tota ! Numerator term of term B for d/h (Eq. 15 Massman et al.)
        real(rk)                   :: fafraczInt_totb ! Denominator term of term B for d/h (Eq. 15 Massman et al.)
        real(rk)                   :: cstress         ! Surface stress at/above canopy height (nondimensional)
        real(rk)                   :: drag            ! Drag area index (i.e., wind speed attenuation) (nondimensional)
        real(rk)                   :: nrat            ! Ratio of drag/cstress (nondimensional)
        real(rk)                   :: z0_set          ! set roughness length (m)
        real(rk)                   :: uc              ! initial guess of wind speed at canopy height (m/s) from log-profile
        real(rk)                   :: lambda_rs        ! local values for influence of roughness sublayer (nondimensional)

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
! Specifically, see Eq. 15 (pg 599) and Eq. 12 (pg 597 of Massman et al.)

        ! Calculate zero-plane displacement height, d/h (Eq. 15 in Massman et al. 2017):
        drag    = CDRAG*PAI
        ! Get first estimate of zpd and z0 to perform initial log wind profile and ustar calculation
        ! Holmes JD. Wind Loading of Structures. 3rd ed. Boca Raton, Florida: CRC Press; 2015.
        ! First assign general z0 dependent on vegetation type and then estimate zpd = 0.75*FCH from Holmes
        if (Z0_OPT .eq. 0) then !Use input model set z0 for first estimate
            z0_set = Z0_MOD
        else if (Z0_OPT .eq. 1) then !Use veg-type dependent z0 first estimate
            if (LU_OPT .eq. 0 .or. LU_OPT .eq. 1) then !VIIRS/MODIS LU types
                !approx/average vegtype mapping to Massman et al. forest types
                if (VTYPE .ge. 1 .and. VTYPE .le. 2) then !VIIRS/MODIS Cat 1-2/Evergreen Needleleaf & Broadleaf
                    z0_set  = 1.0_rk
                end if
                if (VTYPE .ge. 3 .and. VTYPE .le. 5) then !VIIRS/MODIS Cat 3-5/Deciduous Needleleaf, Broadleaf, Mixed Forests
                    z0_set = 1.0_rk
                end if
                if ((VTYPE .ge. 6 .and. VTYPE .le. 10) .or. VTYPE .eq. 12 ) then !VIIRS/MODIS Cat 6-10 or 12/Shrubs, Croplands, and Grasses
                    z0_set = 0.1_rk
                end if
            else
                write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS/MODIS available right now...exiting'
                call exit(2)
            end if
        else
            write(*,*)  'Wrong Z0_OPT choice of ', Z0_OPT, 'in namelist, only 0 or 1 now...exiting'
            call exit(2)
        end if

        !Get wind speed at canopy top with first z0 and zpd=3/4*hc estimates:
        if (HREF > z0_set) then ! input wind speed reference height is > roughness length
            uc = UBZREF*log((FCH-(0.75_rk*FCH)+z0_set)/z0_set)/log(HREF/z0_set)  !MOST
        else                 ! reference height is <= roughness length--at canopy top
            uc = UBZREF
        end if

        !Solve for final z0 and zpd estimates:
        ustrmod = uc*(0.38 - (0.38 + (vonk/log(Z0GHC)))*exp(-1.0*(15.0*drag)))
        cstress = (2.0*(ustrmod**2.0))/(uc**2.0)
        nrat   =  drag/cstress
        qc = 0.60
        qb = 2.0/(1.0 - exp(-1.0))
        qa = 4.04_rk*(-1.0*qb)
        qstar = qa + qb*exp(-1.0*qc*nrat)
        dha =  1.0 - (cosh(qstar*nrat*FCLAI(1))/cosh(qstar*nrat))
        fafraczInt_tota = IntegrateTrapezoid( ZHC,(cosh(qstar*nrat*FCLAI)*ZHC) )
        fafraczInt_totb = IntegrateTrapezoid( ZHC, cosh(qstar*nrat*FCLAI) )
        dhb = fafraczInt_tota/fafraczInt_totb

        ! Final zero-plane displacement (zpd) height
        d_h = dha * dhb

        lambda_rs = LAMBDARS   !set to lambda_rs to user RSL influence term (namelist input LAMBDARS)

        ! Final surface (soil+veg) roughness length, zo/h (Eq. 16 in Massman et al. 2017):
        zo_h  = lambda_rs * (1.0 - d_h) * exp (-vonk*sqrt(2.0/cstress))

    END SUBROUTINE CANOPY_ZPD

end module canopy_profile_mod
