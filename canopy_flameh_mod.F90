module canopy_flameh_mod

    implicit none

contains

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

end module canopy_flameh_mod
