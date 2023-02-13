MODULE canopy_canvars_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Variables Descriptions
! Purpose:  Contains canopy variables descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk

    IMPLICIT NONE

!-------------------------------------------------------------------------------
! Canopy scalars in the model
!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------
! Allocatable canopy variable arrays
!-------------------------------------------------------------------------------

    real(rk), allocatable :: zk            ( : )          ! in-canopy heights (m)
    real(rk), allocatable :: zhc           ( : )          ! z/h
    real(rk), allocatable :: fainc         ( : )          ! incremental foliage shape function
    real(rk), allocatable :: fafracz       ( : )          ! incremental fractional foliage shape function
    real(rk), allocatable :: fafraczInt    ( : )          ! integral of incremental fractional foliage shape function
    real(rk), allocatable :: canBOT        ( : )          ! Canopy bottom wind reduction factors
    real(rk), allocatable :: canTOP        ( : )          ! Canopy top wind reduction factors
    real(rk), allocatable :: canWIND       ( : , : )      ! canopy wind speeds (m/s)
    real(rk), allocatable :: canWIND_3d    ( : , : , : )  ! canopy wind speeds -- 3D (m/s)
    real(rk), allocatable :: dx            ( : )          ! Model grid cell distance/resolution (m)
    real(rk), allocatable :: dx_2d         ( : , : )      ! Model grid cell distance/resolution -- 2D (m)
    real(rk), allocatable :: waf           ( : )          ! Calculated Wind Adjustment Factor
    real(rk), allocatable :: waf_2d        ( : , : )      ! Calculated Wind Adjustment Factor -- 2D
    real(rk), allocatable :: Kz            ( :, : )       ! Eddy Diffusivities (m2/s)
    real(rk), allocatable :: Kz_3d         ( : , : , : )  ! Eddy Diffusivities -- 3D (m2/s)
    real(rk), allocatable :: rjcf          ( :, : )       ! Photolysis Attenuation Correction Factors
    real(rk), allocatable :: rjcf_3d       ( : , : , : )  ! Photolysis Attenuation Correction Factors -- 3D
    real(rk), allocatable :: flameh        ( : )          ! Flame Height (m)
    real(rk), allocatable :: flameh_2d     ( : , : )      ! Flame Height -- 2D (m)

!-------------------------------------------------------------------------------
! Canopy-App Program and version descriptors.
!-------------------------------------------------------------------------------

    CHARACTER(LEN=16),  PARAMETER     :: progname   = 'Canopy-App'
    CHARACTER(LEN=10),  PARAMETER     :: vdate      = '12/22/2022'
    CHARACTER(LEN=8),   PARAMETER     :: ver        = 'V1.0.0'

!-------------------------------------------------------------------------------
! Define output NETCDF data structures.
!-------------------------------------------------------------------------------

    TYPE fld2ddata
        REAL(rk),        POINTER   :: fld        ( : , : )
        CHARACTER(LEN=16)          :: fldname
        CHARACTER(LEN=80)          :: long_name
        CHARACTER(LEN=16)          :: units
        REAL                       :: fillvalue
        CHARACTER(LEN=16)          :: dimnames   ( 4 )
        INTEGER                    :: istart     ( 4 )
        INTEGER                    :: iend       ( 4 )
    END TYPE fld2ddata

    TYPE fld3ddata
        REAL(rk),        POINTER   :: fld        ( : , : , : )
        CHARACTER(LEN=16)          :: fldname
        CHARACTER(LEN=80)          :: long_name
        CHARACTER(LEN=16)          :: units
        REAL                       :: fillvalue
        CHARACTER(LEN=16)          :: dimnames   ( 4 )
        INTEGER                    :: istart     ( 4 )
        INTEGER                    :: iend       ( 4 )
    END TYPE fld3ddata


!-------------------------------------------------------------------------------
! Assign number of time independent and varying 2D/3D fields at cell centers.
!-------------------------------------------------------------------------------

    INTEGER           :: nfld2dxy       ! time-independent 2d cell centers
    INTEGER           :: nfld2dxyt      ! time-varying 2d cell centers
    INTEGER           :: nfld3dxyzt     ! time-varying 3d cell centers

!-------------------------------------------------------------------------------
! Time-independent 2d fields at cell centers.
!-------------------------------------------------------------------------------

    TYPE(fld2ddata), ALLOCATABLE, TARGET :: fld2dxy ( : )
    TYPE(fld2ddata), POINTER     :: g_lat
    TYPE(fld2ddata), POINTER     :: g_lon

!-------------------------------------------------------------------------------
! Time-varying 2d fields at cell centers for output NETCDF
!-------------------------------------------------------------------------------

    TYPE(fld2ddata), ALLOCATABLE, TARGET :: fld2dxyt ( : )
    TYPE(fld2ddata), POINTER     :: c_waf
    TYPE(fld2ddata), POINTER     :: c_flameh

!-------------------------------------------------------------------------------
! Time-varying 3d fields at cell centers for output NETCDF
!-------------------------------------------------------------------------------

    TYPE(fld3ddata), ALLOCATABLE, TARGET :: fld3dxyzt ( : )
    TYPE(fld3ddata), POINTER     :: c_canwind
    TYPE(fld3ddata), POINTER     :: c_Kz
    TYPE(fld3ddata), POINTER     :: c_rjcf

END MODULE canopy_canvars_mod
