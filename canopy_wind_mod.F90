module canopy_wind_mod

    implicit none

contains

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
            call exit(2)
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
            call exit(2)
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

end module canopy_wind_mod
