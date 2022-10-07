module canopy_driver_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_CALCDX(DXOPT, DXSET, NLAT, NLON, LAT, LON, DX )

!-----------------------------------------------------------------------

! Description:
!     computes great circle distance or the orthodromic distance using Haversine
!     formula

! Preconditions:
!     user dx_opt, dx_set, nlat, nlon, lat, and lon

! Subroutines and Functions Called:

! Revision History:
!     Prototype 10/22 by PCC
!     Oct 2022 P.C. Campbell: Initial version
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk       !constants for canopy models
        use canopy_utils_mod, ONLY: CalcDX   !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,    INTENT( IN )  :: DXOPT           ! User DX calculation option
        REAL(RK),   INTENT( IN )  :: DXSET           ! User DX set value if cannot calculate (m)
        INTEGER,    INTENT( IN )  :: NLAT            ! Number of latitude grid cells/points
        INTEGER,    INTENT( IN )  :: NLON            ! Number of longitude grid cells/points
        REAL(RK),   INTENT( IN )  :: LAT(:)          ! Model Latitudes (degrees)
        REAL(RK),   INTENT( IN )  :: LON(:)          ! Model Longitudes (degrees)
        REAL(RK),   INTENT( OUT ) :: DX(:)           ! Distance between two points (m)

!     Local variables
        integer  ::    loc                           ! NLAT*NLON

        do loc=1, NLAT*NLON

            if (DXOPT .eq. 0) then !user set to calculate dx grid cell distance from grid lons
                if (NLON .gt. 1 .and. NLON .gt. 1) then !convert grid points to distances using Haversine formula (m)
                    if (loc .lt. NLAT*NLON) then !inside domain
                        DX(loc) = CalcDX(LAT(loc),LAT(loc+1),LON(loc),LON(loc+1)) !Haversine Calc
                    else !at the domain edge --set to loc-1
                        DX(loc) = DX(NLAT*NLON-1)
                    end if
                else                  !single grid cell/point, use namelist defined dx resolution (m) for cell
                    write(*,*)  'DX_OPT set to calc, but nlon or nlat  <= 1...setting dx = ', &
                        DXSET, ' from namelist'
                    DX(loc) = DXSET
                end if
            else ! user set dx_set from namelist
                DX(loc) = DXSET
            end if

        end do

    END SUBROUTINE CANOPY_CALCDX

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

        !Local variables

!need vegtype to Massman forest mapping

!initializing
        FIRETYPE=-1
        CDRAG=0.0_rk
        PAI=0.0_rk
        ZCANMAX=0.0_rk
        SIGMAU=0.0_rk
        SIGMA1=0.0_rk

        if (LU_OPT .eq. 0) then !VIIRS LU types

            !approx/average vegtype mapping to Massman et al. forest types
            if (VTYPE .ge. 1 .and. VTYPE .le. 2) then !VIIRS Cat 1-2/Evergreen Needleleaf & Broadleaf
                !--> Use average Massman Aspen+Spruce+Pine Forest
                FIRETYPE=0
                CDRAG=(0.20_rk + 0.25_rk + 0.20_rk + 0.20_rk + 0.20_rk)/5.0_rk
                if (PAI_OPT .eq. 0) then      !Katul et al. 2004 vegtype
                    PAI=(5.73_rk + 3.28_rk + 2.41_rk + 2.14_rk + 3.78_rk)/5.0_rk
                else if (PAI_OPT .eq. 1) then !PAI calculation (Massman et al., Eq. 19)
                    PAI=CalcPAI(FCH,FFRAC)
                else if (PAI_OPT .eq. 2) then !PAI = LAI + SAI (WAI)
                    PAI=LAI + 0.52_rk  !WAI  = 0.52 from Toda and Richardson (2018):
                    ! https://doi.org/10.1016/j.agrformet.2017.09.004
                    ! Section 3.3
                else if (PAI_OPT .eq. 3) then !PAI value from user
                    PAI=PAI_SET
                else
                    write(*,*)  'Wrong PAI_OPT choice of ', PAI_OPT, 'in namelist...exiting'
                    call exit(2)
                end if
                ZCANMAX=(0.60_rk + 0.36_rk + 0.60_rk + 0.58_rk + 0.60_rk)/5.0_rk
                SIGMAU=(0.38_rk + 0.60_rk + 0.30_rk + 0.20_rk + 0.10_rk)/5.0_rk
                SIGMA1=(0.16_rk + 0.20_rk + 0.10_rk + 0.20_rk + 0.27_rk)/5.0_rk
            end if

            if (VTYPE .ge. 3 .and. VTYPE .le. 5) then !VIIRS Cat 3-5/Deciduous Needleleaf, Broadleaf, Mixed Forests
                !--> Use Massman Hardwood Forest
                FIRETYPE=0
                CDRAG=0.15_rk
                if (PAI_OPT .eq. 0) then      !Katul et al. 2004 vegtype
                    PAI=4.93_rk
                else if (PAI_OPT .eq. 1) then !Massman PAI calculation (Eq. 19)
                    PAI=CalcPAI(FCH,FFRAC)
                else if (PAI_OPT .eq. 2) then !need PAI function of model LAI
                    PAI=LAI + 0.52_rk  !WAI  = 0.52 from Toda and Richardson (2018):
                    !https://doi.org/10.1016/j.agrformet.2017.09.004
                    ! Section 3.3
                else if (PAI_OPT .eq. 3) then !PAI value from user
                    PAI=PAI_SET
                else
                    write(*,*)  'Wrong PAI_OPT choice of ', PAI_OPT, 'in namelist...exiting'
                    call exit(2)
                end if
                ZCANMAX=0.84_rk
                SIGMAU=0.13_rk
                SIGMA1=0.30_rk
            end if

            if ((VTYPE .ge. 6 .and. VTYPE .le. 10) .or. VTYPE .eq. 12 ) then !VIIRS Cat 6-10 or 12/Shrubs, Croplands, and Grasses
                !--> Average of Massman Corn + Rice )
                FIRETYPE=1
                CDRAG=(0.30_rk + 0.30_rk)/2.0_rk
                if (PAI_OPT .eq. 0) then      !Katul et al. 2004 vegtype
                    PAI=(2.94_rk + 3.10_rk)/2.0_rk
                else if (PAI_OPT .eq. 1) then !PAI calculation (Massman et al., Eq. 19)
                    PAI=CalcPAI(FCH,FFRAC)
                else if (PAI_OPT .eq. 2) then !PAI = LAI + SAI (WAI)
                    PAI=LAI + 0.52_rk  !WAI  = 0.52 from Toda and Richardson (2018):
                    !https://doi.org/10.1016/j.agrformet.2017.09.004
                    ! Section 3.3
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
            write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS available right now...exiting'
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
        fatot             = 0.0_rk
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

        if(allocated(fainc))      deallocate(fainc)
        if(allocated(fafracz))    deallocate(fafracz)

    END SUBROUTINE CANOPY_FOLIAGE

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_ZPD( ZHC, FCLAI, UBZREF, Z0GHC, &
        LAMDARS, CDRAG, PAI, d_h, zo_h )

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
        REAL(RK),    INTENT( IN )  :: UBZREF          ! Mean wind speed at zref-height of canopy top (m/s)
        REAL(RK),    INTENT( IN )  :: Z0GHC           ! Ratio of ground roughness length to canopy top height (nondimensional)
        REAL(RK),    INTENT( IN )  :: LAMDARS         ! Value representing influence of roughness sublayer (nondimensional)
        REAL(RK),    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
        REAL(RK),    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
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
        real(rk)                   :: cstress         ! Suface stress at/above canopy height (nondimensional)
        real(rk)                   :: drag            ! Drag area index (i.e., wind speed attentuation) (nondimensional)
        real(rk)                   :: nrat            ! Ratio of drag/cstress (nondimensional)

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
! Specifically, see Eq. 15 (pg 599) and Eq. 12 (pg 597 of Massman et al.)

        ! Calculate zero-plane displacement height, d/h (Eq. 15 in Massman et al. 2017):
        drag    = CDRAG*PAI
        ustrmod = UBZREF*(0.38 - (0.38 + (vonk/log(Z0GHC)))*exp(-1.0*(15.0*drag)))
        cstress = (2.0*(ustrmod**2.0))/(UBZREF**2.0)
        nrat   =  drag/cstress
        qc = 0.60
        qb = 2.0/(1.0 - exp(-1.0))
        qa = 4.04*(-1.0*qb)
        qstar = qa + qb*exp(-1.0*qc*nrat)
        dha =  1.0 - (cosh(qstar*nrat*FCLAI(1))/cosh(qstar*nrat))
        fafraczInt_tota = IntegrateTrapezoid( ZHC,(cosh(qstar*nrat*FCLAI)*ZHC) )
        fafraczInt_totb = IntegrateTrapezoid( ZHC, cosh(qstar*nrat*FCLAI) )
        dhb = fafraczInt_tota/fafraczInt_totb

        ! zero-plane displacement height
        d_h = dha * dhb
        ! Calculate surface (soil+veg) roughness length, zo/h (Eq. 16 in Massman et al. 2017):
        zo_h  = LAMDARS * (1.0 - d_h) * exp (-vonk*sqrt(2.0/cstress))

    END SUBROUTINE CANOPY_ZPD

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_WIND( HCM, ZK, FAFRACK, UBZREF, Z0GHC, &
        CDRAG, PAI, HREF, D_H, ZO_H, MOL, RSL_OPT, &
        CANBOT_OUT, CANTOP_OUT, CANWIND )

!-----------------------------------------------------------------------

! Description:
!     computes mean wind speed for given height (z) below the canopy top.

! Preconditions:
!     in-canopy height, at canopy wind speed, plant distribution functions

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone canopy wind model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk, vonk  !constants for canopy models

! Arguments:
!       IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: ZK              ! Above/Below canopy height, z (m)
        REAL(RK),    INTENT( IN )  :: FAFRACK         ! Fractional (z) shapes of the
        ! plant surface distribution (nondimensional)
        REAL(RK),    INTENT( IN )  :: UBZREF          ! Mean wind speed at reference height (m/s)
        REAL(RK),    INTENT( IN )  :: Z0GHC          ! Ratio of ground roughness length to canopy top height (nondimensional)
        REAL(RK),    INTENT( IN )  :: CDRAG           ! Drag coefficient (nondimensional)
        REAL(RK),    INTENT( IN )  :: PAI             ! Total plant/foliage area index (nondimensional)
        REAL(RK),    INTENT( IN )  :: HREF            ! Reference Height above canopy asssociated with ref wind speed  (m)
        REAL(RK),    INTENT( IN )  :: D_H             ! Zero plane displacement heights (nondimensional)
        REAL(RK),    INTENT( IN )  :: ZO_H            ! Surface (soil+veg) roughness lengths (nondimensional)
        REAL(RK),    INTENT( IN )  :: MOL             ! Model input Monin-Obukhov Length
        INTEGER,     INTENT( IN )  :: RSL_OPT         ! RSL option used in model from Rosenzweig et al. 2021 (default = 0, off)
        REAL(RK),    INTENT( OUT ) :: CANBOT_OUT      ! Canopy bottom wind reduction factor = canbot (nondimensional)
        REAL(RK),    INTENT( OUT ) :: CANTOP_OUT      ! Canopy top wind reduction factor = cantop    (nondimensional)
        REAL(RK),    INTENT( OUT ) :: CANWIND         ! Mean canopy wind speed at current z (m/s)
!       Local variables
        real(rk)                   :: ustrmod         ! Friction Velocity parameterization based on Massman 2017 (m/s)
        real(rk)                   :: z0g             ! Ground roughness length based on z0g/HCCM ratio (m)
        real(rk)                   :: zkhcm           ! Current zk/hcm ratio (nondimensional)
        real(rk)                   :: cstress         ! Suface stress at/above canopy height (nondimensional)
        real(rk)                   :: drag            ! Drag area index (i.e., wind speed attentuation) (nondimensional)
        real(rk)                   :: nrat            ! Ratio of drag/cstress (nondimensional)
        real(rk)                   :: canbot          ! Logarithmic wind speed that is dominant near the ground (nondimensional)
        real(rk)                   :: cantop          ! Hyperbolic cosine wind speed that is dominant near the top of canopy (nondimensional)
        real(rk)                   :: zpd             ! Zero plane displacement heights (m)
        real(rk)                   :: z0m             ! Surface (soil+veg) roughness lengths (m)
        real(rk)                   :: uc              ! Wind directly at canopy top (m/s)
        real(rk)                   :: beta, beta_n    ! Ratio of u* to wind speed at canopy height
        real(rk)                   :: beta1, beta2    ! Terms used in calculating beta above
        real(rk)                   :: can_length      ! Canopy length scale (= HCM/CDRAG*PAI)

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354
        zkhcm = ZK/HCM
        z0g = Z0GHC*HCM

        ! Nondimensional canopy wind speed term that dominates near the ground:
        if (ZK >= z0g .and. ZK <= HCM) then
            canbot = log(zkhcm/Z0GHC)/log(1.0_rk/Z0GHC)
        else if (ZK >= 0 .and. ZK <= z0g) then
            canbot = 0.0  !No-slip condition at surface (u=0 at z=0)
        else
            canbot = 1  ! should be zk > hcm
        end if

        CANBOT_OUT = canbot
        ! Nondimensional canopy wind speed term that dominates near the top of the canopy:
        ! Assume the drag area distribution over depth of canopy can be approx. p1=0 (no shelter factor) and d1=0
        ! (no drag coefficient relation to wind speed) -- thus no intergration then required in Eq. (4) of Massman et al.
        drag    = CDRAG*PAI

        !first adjust reference wind down to canopy top wind using MOST (with no RSL effects)
        zpd = D_H*HCM  !zero-plane displacement height (not scaled to HCM)
        z0m = ZO_H*HCM !aerodynamic roughness lengh (not scaled to HCM)
        if((HCM-zpd) <= 0.) then
            write(*,*) "critical problem: hcan <= zpd"
            call exit(1)
        end if

        if (HREF > z0m) then ! input wind speed reference height is > roughness length
            uc = UBZREF*log((HCM-zpd+z0m)/z0m)/log(HREF/z0m)  !MOST
        else                 ! reference height is <= roughness length--at canopy top
            uc = UBZREF
        end if

        !calculate U* from Massman 1997 (https://doi.org/10.1023/A:1000234813011)
        ustrmod = uc*(0.38_rk - (0.38_rk + (vonk/log(Z0GHC)))*exp(-1.0_rk*(15.0_rk*drag)))

        if (HREF > z0m) then ! input wind speed reference height is > roughness length
            if (RSL_OPT .eq. 0) then !MOST From NoahMP (M. Barlarge)
                uc = uc  !no RSL effects--Just MOST
            else if (RSL_OPT .eq. 1) then !Unified RSL (Rosenzweig et al., 2021)  https://doi.org/10.1029/2021MS002665
                beta_n = 0.35_rk               !beta for neutral conditions
                can_length = HCM/drag          !canopy length scale

                !Check stability for β calculation
                IF ( MOL .LT.  0.0 )  THEN     !UNSTABLE
                    beta1 = -1.0_rk*((16.0_rk*(beta_n**4.0_rk)*can_length)/MOL)
                    beta2 = sqrt( (((16.0_rk*(beta_n**4.0_rk)*can_length)/MOL)**2.0_rk) + (4.0_rk*(beta_n**4.0_rk)) )
                    beta = sqrt( 0.5_rk * (beta1 + beta2) )
                ELSE IF ( MOL .GT.  0.0 ) THEN !STABLE to VERY STABLE
                    beta1 = -1.0_rk * (( 0.5 * ( sqrt ((-27.0_rk*beta_n*((5.0_rk*can_length)/MOL)**2.0_rk)**2.0_rk + &
                        4.0_rk*(15.0_rk*(can_length/MOL))**3.0_rk ) - &
                        27.0_rk*beta_n*((5.0_rk*can_length)/MOL)**2.0_rk )) ** (1.0_rk/3.0_rk)) * &
                        (MOL/(15.0_rk*can_length))
                    beta2 = +1.0_rk * (( 0.5 * ( sqrt ((-27.0_rk*beta_n*((5.0_rk*can_length)/MOL)**2.0_rk)**2.0_rk + &
                        4.0_rk*(15.0_rk*(can_length/MOL))**3.0_rk ) - &
                        27.0_rk*beta_n*((5.0_rk*can_length)/MOL)**2.0_rk )) ** (-1.0_rk/3.0_rk))
                    beta = (beta1 + beta2)
                ELSE                          !NEUTRAL
                    beta = beta_n
                END IF

                !The solution for β is constrained between 0.2 and 0.5 (Bonan et al., 2018).
                if (beta .lt. 0.2) then
                    beta = 0.2
                else if (beta .gt. 0.5) then
                    beta = 0.5
                end if

                uc = ustrmod/beta
            else
                write(*,*) 'wrong namelist option', RSL_OPT, 'only 0 (MOST) or 1 (RSL) available'
                call exit(2)
            end if
        else
            uc = UBZREF
        end if

        if (uc > UBZREF) then !reference height still small and close to roughness length
            uc = UBZREF
        end if

        if (uc <= 0.0) then
            write(*,*)  'Uc cannot be <= 0 at  ',uc , ' in canopy_wind calc...exiting'
            call exit(1)
        end if

        !Calculate remaining in-canopy parameters (Massman et al., 2017)
        cstress = (2.0_rk*(ustrmod**2.0_rk))/(uc**2.0_rk)
        nrat   =  drag/cstress
        cantop = cosh(nrat*FAFRACK)/cosh(nrat)
        CANTOP_out = cantop

        if (ZK <= HCM) then       !at or below canopy top --> modify uc winds
            CANWIND = uc*canbot*cantop
        else
            CANWIND = UBZREF       !above canopy top constant and equal to reference winds
        end if

    END SUBROUTINE CANOPY_WIND

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_FLAMEH( FLAMEH_OPT, FLAMEH_SET, DX, MODRES, &
        FRP, MIDFLAMEPOINT, FLAMEH )

!-----------------------------------------------------------------------

! Description:
!     computes flame height and midflamepoint layer

! Preconditions:
!     user flameh option and set, grid resolution dx, vert model res, and FRP

! Subroutines and Functions Called:

! Revision History:
!     Oct 2022 P.C. Campbell: Initial flame height calculation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk         !constants for canopy models
        use canopy_utils_mod, ONLY: CalcFlameH !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )  :: FLAMEH_OPT      ! Integer for flameh values used or calculated (default = 0)
        REAL(RK),    INTENT( IN )  :: FLAMEH_SET      ! User Set Flame Height (m)
        REAL(RK),    INTENT( IN )  :: DX              ! DX cell distances using haversine formula (m)
        REAL(RK),    INTENT( IN )  :: MODRES          ! Canopy model input vertical resolution (m)
        REAL(RK),    INTENT( IN )  :: FRP             ! Model input FRP (MW/grid cell area)
        INTEGER,     INTENT( OUT ) :: MIDFLAMEPOINT   ! Indice of the mid-flame point
        REAL(RK),    INTENT( OUT ) :: FLAMEH          ! Flame Height (m)

!     Local variables

        integer  ::    flamelays                      ! number of flame layers

        if (FLAMEH_OPT .eq. 0) then
            if (DX .le. 0.0) then !single lon point -- reset to user value
                write(*,*)  'FLAMEH_OPT set to FRP calc, but dx <= 0...setting flameh = ', &
                    FLAMEH_SET, ' from namelist'
                FLAMEH = FLAMEH_SET
            else                                    !calculate flameh
                FLAMEH = CalcFlameH(FRP,DX)
            end if
        else if (FLAMEH_OPT .eq. 1) then  !user set value
            FLAMEH = FLAMEH_SET
        else
            write(*,*)  'Wrong FLAMEH_OPT choice of ', FLAMEH_OPT, ' in namelist...exiting'
            call exit(2)
        end if
        if (FLAMEH .lt. MODRES) then !flameh beween first (z=0) and second layer height
            MIDFLAMEPOINT = 2    !do not put at z = 0, but rather in second layer
        else
            flamelays     = floor(FLAMEH/MODRES)
            MIDFLAMEPOINT = max((flamelays/2),2)
        end if

    END SUBROUTINE CANOPY_FLAMEH

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_WAF( HCM, LAMDARS, HREF, FLAMEH, FIRETYPE, &
        D_H, ZO_H, CANBOTMID, CANTOPMID, WAF )

!-----------------------------------------------------------------------

! Description:
!     computes Wind Adjustment Factor for fire spread for either sub- or above-canopy fires.

! Preconditions:
!     in-canopy height, firetype and height, mean mid-canopy wind speed, and plant distribution functions

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Massman et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone canopy wind model
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk !constants for canopy models
        use canopy_utils_mod           !utilities for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: LAMDARS         ! Value representing influence of roughness sublayer (nondimensional)
        REAL(RK),    INTENT( IN )  :: HREF            ! Reference Height (m) above the canopy
        REAL(RK),    INTENT( IN )  :: FLAMEH          ! Flame Height (m) -- Only for Above Canopy Fire
        INTEGER ,    INTENT( IN )  :: FIRETYPE        ! 1 = Above Canopy Fire; 0 = Below Canopy Fire
        REAL(RK),    INTENT( IN )  :: CANBOTMID       ! Mid-flame canopy bottom wind reduction factor (nondimensional)
        REAL(RK),    INTENT( IN )  :: CANTOPMID       ! Mid-flame canopy top wind reduction factor (nondimensional)
        REAL(RK),    INTENT( IN )  :: D_H             ! Zero-plane displacement height, d/h
        REAL(RK),    INTENT( IN )  :: ZO_H            ! Surface (soil+veg) roughness length, zo/h
        REAL(RK),    INTENT( OUT ) :: WAF             ! Wind Adjustment Factor (nondimensional)
!     Local variables
        real(rk)                   :: term1           ! Major Term1 in WAF calculation (Eqs. 17 and 18 Massman et al. 2017)
        real(rk)                   :: term2           ! Major Term2 in WAF calculation (Eqs. 17 and 18 Massman et al. 2017)
        real(rk)                   :: delta           ! Ratio parmeter used in WAF for above-canopy (Eq. 18 Massman et al.)

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354

        if (HREF <= 0.0) then
            write(*,*) "critical problem: HREF <= 0, WAF calculation not accurate and thus is reset to 1"
            waf = 1.0
        else
            if (FIRETYPE == 0) then  !sub-canopy
                term1 = log( LAMDARS * ( (1.0 - D_H)/ZO_H ) )  !numerator
                term2 = log( LAMDARS * ( ( (HREF/HCM) + 1.0 - D_H ) / ZO_H ) )  !denominatory
                waf   = CANBOTMID * CANTOPMID * (term1 / term2)
            else                     !above-canopy
                delta = (1.0 - D_H) / (FLAMEH/HCM)
                term1 = log( LAMDARS * ( ( (FLAMEH/HCM) +  1.0 - D_H ) / ZO_H ) ) - &   !numerator
                    1.0 + (delta*log((1.0/delta) + 1.0))
                term2 = log( LAMDARS * ( ( ( (HREF/HCM) + 1.0 - D_H) )/ ZO_H ) )  !denominator
                waf   = term1 / term2
            end if
        end if

    END SUBROUTINE CANOPY_WAF

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_PHOT( FCLAI, LAI, CLU, COSZEN, RJCF )

!-----------------------------------------------------------------------

! Description:
!     computes photolysis attenuation in canopy.

! Preconditions:
!     in-canopy height, and model LAI, clumping index, and solar zenith angle

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Makar et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone vertical diffusion calculation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk     !constants for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: FCLAI(:)           ! Model input Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )  :: LAI             ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )  :: CLU             ! Model input Clumping Index
        REAL(RK),    INTENT( IN )  :: COSZEN          ! Model input Cosine Solar Zenith Angle
        REAL(RK),    INTENT( OUT ) :: RJCF(:)          ! Photolysis correction factor

!     Local variables

! Citation:
!Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone.
!Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243

!Eq. 1 on pg 9 of Makar et al., 2017
        RJCF = MAX(1.0E-10_rk, EXP(-1.0_rk*(0.5_rk*(LAI*(1.0_rk-FCLAI))*CLU)/MAX(0.05_rk, COSZEN)))

    END SUBROUTINE CANOPY_PHOT
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE CANOPY_EDDYX( HCM, ZK, USTAR, MOL, KZ )

!-----------------------------------------------------------------------

! Description:
!     computes Eddy Diffusivity, Kz, within and just above canopy.

! Preconditions:
!     in-canopy height, and model friction velocity and MOL

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/22 by PCC, based on Makar et al. (2017) algorithms
!     Jun 2022 P.C. Campbell: Initial standalone vertical diffusion calculation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk, pi !constants for canopy models

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )  :: HCM             ! Height of canopy top (m)
        REAL(RK),    INTENT( IN )  :: ZK              ! Above/Below canopy height, z (m)
        REAL(RK),    INTENT( IN )  :: USTAR           ! Model input friction velocity (m/s)
        REAL(RK),    INTENT( IN )  :: MOL             ! Model input Monin-Obukhov Length
        REAL(RK),    INTENT( OUT ) :: KZ              ! Estimated Eddy Diffusivity with canopy turbulence

!     Local variables
        real(rk)                   :: hol             ! local canopy stability parameter (hc/MOL)
        real(rk)                   :: tlc             ! turbulence length scale in canopy (m)
        real(rk)                   :: rr              ! sigma parameter R (Makar et al., 2017)
        real(rk)                   :: aa              ! sigma parameter A (Makar et al., 2017)
        real(rk)                   :: bb              ! sigma parameter B (Makar et al., 2017)
        real(rk)                   :: sigma           ! eulerian vertial velocity (m/s)

! Citation:
!Makar, P., Staebler, R., Akingunola, A. et al. The effects of forest canopy shading and turbulence on boundary layer ozone.
!Nat Commun 8, 15243 (2017). https://doi.org/10.1038/ncomms15243
!Raupauch M. R. A Practical Lagrangian method for relating scalar concentrations to
! source distributions in vegetation canopies. Q. J. R. Meteor. Soc. (1989), 115, pp 609-632
! Eqs. 2-9 on pgs 10-11.

        hol = HCM/MOL
        tlc = (HCM/USTAR) * (                        &
            (0.256_rk * (ZK-(0.75_rk*HCM))/HCM ) +      &
            (0.492_rk*EXP((-0.256_rk*ZK/HCM)/0.492_rk)) )
        sigma = 0.0_rk

        IF ( hol .LT. -0.1 )  THEN  !UNSTABLE
            IF ( ZK/HCM .GT. 1.25_rk ) THEN !SIGMACAN = Eulerian vertical velocity variance
                sigma = 1.25_rk*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                sigma = USTAR * ( 0.75_rk + (0.5_rk * COS((pi/1.06818_rk) *     &
                    (1.25_rk - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25_rk*USTAR
            END IF
        END IF

        IF ( hol .GE. -0.1  .AND. hol .LT. 0.1 )   THEN  !NEUTRAL
            IF ( ZK/HCM .GT. 1.25 ) THEN
                sigma = 1.0_rk*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                sigma = USTAR * ( 0.625_rk + (0.375_rk * COS((pi/1.06818_rk) *  &
                    (1.25_rk - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25_rk*USTAR
            END IF
        END IF

        IF ( hol .GE.  0.1  .AND.  hol .LT. 0.9 )   THEN  !STABLE
            IF ( ZK/HCM .GT. 1.25 ) THEN
                sigma = 0.25_rk*(4.375_rk - (3.75_rk*hol))*USTAR
            END IF
            IF ( ZK/HCM .GE. 0.175  .AND.  ZK/HCM .LE. 1.25 ) THEN
                rr=4.375_rk-(3.75_rk*hol)
                aa=(0.125_rk*rr) + 0.125_rk
                bb=(0.125_rk*rr) - 0.125_rk
                sigma = USTAR * ( aa + (bb * COS((pi/1.06818_rk) *   &
                    (1.25_rk - (ZK/HCM)))) )
            END IF
            IF ( ZK/HCM .LT. 0.175 )  THEN
                sigma = 0.25_rk*USTAR
            END IF
        END IF

        IF ( hol .GE.  0.9 ) THEN  !VERY STABLE
            sigma = 0.25_rk*USTAR
        END IF

        KZ = (sigma*sigma)*tlc

    END SUBROUTINE CANOPY_EDDYX

end module canopy_driver_mod
