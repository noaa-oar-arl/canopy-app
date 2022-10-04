MODULE canopy_canmet_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Met/Sfc Input Variable Descriptions
! Purpose:  Contains canopy met and sfc input  variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------

    IMPLICIT NONE

!! .... defines canopy optionss (read from user namelist)
    INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
! Generic 2D met/sfc input variables that should be passed to canopy calculations
    TYPE :: variable_type
        real(rk)   :: lat          !latitude of cell/point
        real(rk)   :: lon          !longitude of cell/point
        real(rk)   :: fh           !forest/canopy height
        real(rk)   :: ws           !wind speed at reference height above canopy (10 m)
        real(rk)   :: clu          !clumping index
        real(rk)   :: lai          !leaf area index
        integer    :: vtype        !vegetation type
        real(rk)   :: ffrac        !forest fraction
        real(rk)   :: ust          !friction velocity (u*)
        real(rk)   :: csz          !cosine of solar zenith angle
        real(rk)   :: z0           !surface roughness length
        real(rk)   :: mol          !Monin-Obukhov length
        real(rk)   :: frp          !fire radiative power (MW per grid cell)
    end TYPE variable_type

    type(variable_type), allocatable :: variables( : )

    ! Met/Sfc variable reassignment names  above reference conditions from the model
    real(rk)       ::    hcmref          !Input Canopy Height (m)
    real(rk)       ::    ubzref          !Input above canopy/reference 10-m model wind speed (m/s)
    real(rk)       ::    cluref          !Input canopy clumping index
    real(rk)       ::    lairef          !Input leaf area index
    integer        ::    vtyperef        !Input vegetation type (VIIRS)
    real(rk)       ::    ffracref        !Input forest fraction of grid cell
    real(rk)       ::    ustref          !Input friction velocity
    real(rk)       ::    cszref          !Input cosine of zenith angle
    real(rk)       ::    z0ref           !Input total/surface roughness length
    real(rk)       ::    molref          !Input Monin-Obukhov Length
    real(rk)       ::    frpref          !Input fire radiative power

END MODULE canopy_canmet_mod
