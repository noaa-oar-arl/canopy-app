MODULE canopy_canvars_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Variables Descriptions
! Purpose:  Contains canopy variables descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------

    IMPLICIT NONE

!! .... defines canopy variables calculated in the model
    INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 307)
    integer        ::    firetype      !1 = Above Canopy Fire; 0 = Below Canopy Fire
    integer        ::    cansublays    !Number of sub-canopy layers
    integer        ::    midflamepoint !Indice of the mid-flame point
    real(rk)       ::    cdrag         !Drag coefficient (nondimensional)
    real(rk)       ::    pai           !Plant/foliage area index (nondimensional)
    real(rk)       ::    zcanmax       !Height of maximum foliage area density (z/h) (nondimensional)
    real(rk)       ::    sigmau        !Standard deviation of shape function above zcanmax (z/h)
    real(rk)       ::    sigma1        !Standard deviation of shape function below zcanmax (z/h)
    real(rk)       ::    d_h           !Zero plane displacement heights (z/h)
    real(rk)       ::    zo_h          !Surface (soil+veg) roughness lengths (z/h)
    real(rk)       ::    flameh        !Flame Height (m)

!allocatable canopy variables
    real(rk), allocatable :: zk         ( : )     ! in-canopy heights (m)
    real(rk), allocatable :: zhc        ( : )     ! z/h
    real(rk), allocatable :: fainc      ( : )     ! incremental foliage shape function
    real(rk), allocatable :: fafracz    ( : )     ! incremental fractional foliage shape function
    real(rk), allocatable :: fafraczInt ( : )     ! integral of incremental fractional foliage shape function
    real(rk), allocatable :: canBOT     ( : )     ! Canopy bottom wind reduction factors
    real(rk), allocatable :: canTOP     ( : )     ! Canopy top wind reduction factors
    real(rk), allocatable :: canWIND    ( :, : )  ! canopy wind speeds (m/s)
    real(rk), allocatable :: dx         ( : )     ! Model grid cell distance/resolution (m)
    real(rk), allocatable :: waf        ( : )     ! Calculated Wind Adjustment Factor
    real(rk), allocatable :: Kz         ( :, : )  ! Eddy Diffusivities (m2/s)
    real(rk), allocatable :: rjcf       ( :, : )  ! Photolysis Attenuation Correction Factors


END MODULE canopy_canvars_mod
