

# File canopy\_fire\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_fire\_mod.F90**](canopy__fire__mod_8F90.md)

[Go to the documentation of this file](canopy__fire__mod_8F90.md)


```Fortran



module canopy_fire_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_flameh( FLAMEH_OPT, FLAMEH_SET, DX, MODRES, &
        FRP_IN, FRP_FAC, FCH, LU_OPT, VTYPE, FLAMEH_CAL, &
        MIDFLAMEPOINT, FLAMEH )

!-----------------------------------------------------------------------

! Description:
!     computes flame height and midflamepoint layer needed for WAF

! Preconditions:
!     user flameh option and set, grid resolution dx, vert model res, and FRP

! Subroutines and Functions Called:

! Revision History:
!     Oct 2022 P.C. Campbell: Initial flame height calculation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk         !constants for canopy models
        use canopy_utils_mod, ONLY: calcflameh !utilities for canopy models

! Arguments:
!     IN/OUT
        INTEGER,     INTENT( IN )  :: FLAMEH_OPT      ! Integer for flameh values used or calculated (default = 0)
        REAL(RK),    INTENT( IN )  :: FLAMEH_SET      ! User Set Flame Height (m)
        REAL(RK),    INTENT( IN )  :: DX              ! DX cell distances using haversine formula (m)
        REAL(RK),    INTENT( IN )  :: MODRES          ! Canopy model input vertical resolution (m)
        REAL(RK),    INTENT( IN )  :: FRP_IN          ! Model input FRP (MW/grid cell area)
        REAL(RK),    INTENT( IN )  :: FRP_FAC         ! FRP tuning factor for flame height calculation
        REAL(RK),    INTENT( IN )  :: FCH             ! Grid cell canopy height (m)
        INTEGER,     INTENT( IN )  :: LU_OPT          ! Supported land use classifications
        INTEGER,     INTENT( IN )  :: VTYPE           ! Dominant vegetation type
        INTEGER,     INTENT( IN )  :: FLAMEH_CAL      ! Option of vegtype dependent FRP to Flame Height relationships used
        INTEGER,     INTENT( OUT ) :: MIDFLAMEPOINT   ! Indice of the mid-flame point
        REAL(RK),    INTENT( OUT ) :: FLAMEH          ! Flame Height (m)

!     Local variables

        integer  ::    flamelays

        real(rk) ::    frp

        frp = frp_in*frp_fac  !apply FRP tuning factor for flame height

        if (flameh_opt .eq. 0) then
            if (dx .le. 0.0) then !single lon point -- reset to user value
                write(*,*)  'FLAMEH_OPT set to FRP calc, but dx <= 0...setting flameh = ', &
                    flameh_set, ' from namelist'
                flameh = flameh_set
            else                                    !calculate flameh
                flameh = calcflameh(frp,dx,lu_opt,vtype,flameh_cal)
            end if
        else if (flameh_opt .eq. 1) then  !user set value
            flameh = flameh_set
        else if (flameh_opt .eq. 2) then  !both FRP calc and user set
            if (frp .gt. 0.0) then
                flameh = calcflameh(frp,dx,lu_opt,vtype,flameh_cal)
            else
                flameh = flameh_set
            end if
        else if (flameh_opt .eq. 3) then  !overides use of "flame heights" --> User set fraction (<= 1) of FCH
            if (flameh_set .le. 1.0) then
                flameh = fch * flameh_set  !not real flame height but uses WAF at this fractional FCH
            else
                write(*,*)  'Under this FLAMEH_OPT ', flameh_opt, ' FLAMEH_SET must be =< 1.0'
                call exit(2)
            end if
        else if (flameh_opt .eq. 4) then  !uses FRP and overide elsewhere
            if (frp .gt. 0.0) then
                flameh = calcflameh(frp,dx,lu_opt,vtype,flameh_cal)
            else
                if (flameh_set .le. 1.0) then
                    flameh = fch * flameh_set  !not real flame height but uses WAF at this fractional FCH
                else
                    write(*,*)  'Under this FLAMEH_OPT ', flameh_opt, ' FLAMEH_SET must be =< 1.0'
                    call exit(2)
                end if
            end if
        else if (flameh_opt .eq. 5) then  !uses adjusted FRP (based on intensity) and overide elsewhere.
            if (frp .gt. 0.0) then
                if ( ((frp*1000.0_rk)/dx) .ge. 1700.0_rk ) then       !Crown fire likely (Andrews et al., 2011).
                    flameh = fch                                      !https://doi.org/10.2737/RMRS-GTR-253
                else
                    flameh = calcflameh(frp,dx,lu_opt,vtype,flameh_cal)
                end if
            else
                if (flameh_set .le. 1.0) then
                    flameh = fch * flameh_set  !not real flame height but uses WAF at this fractional FCH
                else
                    write(*,*)  'Under this FLAMEH_OPT ', flameh_opt, ' FLAMEH_SET must be =< 1.0'
                    call exit(2)
                end if
            end if
        else
            write(*,*)  'Wrong FLAMEH_OPT choice of ', flameh_opt, ' in namelist...exiting'
            call exit(2)
        end if
        if (flameh .le. modres) then !flameh between first (z=0) and second layer height
            midflamepoint = 2    !do not put at z = 0, but rather in second layer
        else
            flamelays     = floor(flameh/modres) + 1      !force full flame layers
            midflamepoint = max(ceiling(flamelays/2.0),2) !conservative midflamepoint
            if ( flameh_opt .eq. 5 .and. ((frp*1000.0_rk)/dx) .ge. 1700.0_rk ) then !crowning likely
                midflamepoint = max(ceiling(flamelays/1.0),2) !sets to top of flame
            end if
        end if

    END SUBROUTINE canopy_flameh

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_waf( HCM, LAMBDARS, HREF, FLAMEH, FIRETYPE, &
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
        REAL(RK),    INTENT( IN )  :: LAMBDARS         ! Value representing influence of roughness sublayer (nondimensional)
        REAL(RK),    INTENT( IN )  :: HREF            ! Reference Height (m) above the canopy
        REAL(RK),    INTENT( IN )  :: FLAMEH          ! Flame Height (m) -- Only for Above Canopy Fire
        INTEGER ,    INTENT( IN )  :: FIRETYPE        ! 1 = Above Canopy Fire; 0 = Below Canopy Fire
        REAL(RK),    INTENT( IN )  :: CANBOTMID       ! Mid-flame canopy bottom wind reduction factor (nondimensional)
        REAL(RK),    INTENT( IN )  :: CANTOPMID       ! Mid-flame canopy top wind reduction factor (nondimensional)
        REAL(RK),    INTENT( IN )  :: D_H             ! Zero-plane displacement height, d/h
        REAL(RK),    INTENT( IN )  :: ZO_H            ! Surface (soil+veg) roughness length, zo/h
        REAL(RK),    INTENT( OUT ) :: WAF             ! Wind Adjustment Factor (nondimensional)

!     Local variables
        real(rk)                   :: term1

        real(rk)                   :: term2

        real(rk)                   :: delta

        real(rk)                   :: lambda_rs

! Citation:
! An improved canopy wind model for predicting wind adjustment factors and wildland fire behavior
! (2017)  W.J. Massman, J.M. Forthofer, M.A. Finney.  https://doi.org/10.1139/cjfr-2016-0354

        if (href <= 0.0) then
            write(*,*) "warning: HREF <= 0, WAF calculation not accurate and thus is reset to 1"
            waf = 1.0
        else

            lambda_rs = lambdars

            if (firetype == 0) then  !sub-canopy
                term1 = log( lambda_rs * ( (1.0 - d_h)/zo_h ) )  !numerator
                term2 = log( lambda_rs * ( ( (href/hcm) + 1.0 - d_h ) / zo_h ) )  !denominatory
                waf   = canbotmid * cantopmid * (term1 / term2)
            else                     !above-canopy
                delta = (1.0 - d_h) / (flameh/hcm)
                term1 = log( lambda_rs * ( ( (flameh/hcm) +  1.0 - d_h ) / zo_h ) ) - &   !numerator
                    1.0 + (delta*log((1.0/delta) + 1.0))
                term2 = log( lambda_rs * ( ( ( (href/hcm) + 1.0 - d_h) )/ zo_h ) )  !denominator
                waf   = term1 / term2
            end if
        end if

    END SUBROUTINE canopy_waf

end module canopy_fire_mod
```


