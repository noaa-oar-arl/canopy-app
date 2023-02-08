MODULE canopy_canmet_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Met/Sfc Input Variable Descriptions
! Purpose:  Contains canopy met and sfc input  variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk

    IMPLICIT NONE
!! .... defines canopy options (read from user namelist)

! Generic 2D met/sfc input variables that should be passed to canopy calculations
    TYPE :: variable_type
        real(rk)   :: lat          !latitude of cell/point
        real(rk)   :: lon          !longitude of cell/point
        real(rk)   :: fh           !forest/canopy height (m)
        real(rk)   :: ugrd10m      !u wind speed at reference height above canopy (m/s)
        real(rk)   :: vgrd10m      !v wind speed at reference height above canopy (m/s)
        real(rk)   :: clu          !clumping index
        real(rk)   :: lai          !leaf area index
        integer    :: vtype        !vegetation type
        real(rk)   :: ffrac        !forest fraction
        real(rk)   :: fricv        !friction velocity (u*) (m/s)
        real(rk)   :: csz          !cosine of solar zenith angle
        real(rk)   :: sfcr         !surface roughness length (m)
        real(rk)   :: mol          !Monin-Obukhov length (m)
        real(rk)   :: frp          !fire radiative power (MW)
        real(rk)   :: href         !reference height above the canopy (m)
    end TYPE variable_type

    type(variable_type), allocatable :: variables( : ), variables_2d( : , :)

    ! Met/Sfc variable reassignment names  above reference conditions from the model
    real(rk)       ::    latref          !latitude of cell/point
    real(rk)       ::    lonref          !longitude of cell/point
    real(rk)       ::    hcmref          !Input Canopy Height
    real(rk)       ::    uref            !Input above canopy/reference 10-m U wind speed
    real(rk)       ::    vref            !Input above canopy/reference 10-m V wind speed
    real(rk)       ::    ubzref          !Input above canopy/reference 10-m model wind speed
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
