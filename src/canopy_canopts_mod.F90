MODULE canopy_canopts_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Option Variable Descriptions
! Purpose:  Contains canopy option variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk
    IMPLICIT NONE

!! .... defines canopy options (read from user namelist)
    integer             ::    infmt_opt   !Integer for choosing 1D or 2D input file format (default = 0, 2D)
    integer             ::    href_opt    !Integer for using set href in namelist (=0) or array from file(=1) (default = 0)
    real(rk)            ::    href_set    !Set reference Height above canopy @ 10 m  (m)
    logical             ::    ifcanwind   !logical canopy wind option (default = .FALSE.)
    logical             ::    ifcanwaf    !logical canopy WAF option (default = .FALSE.)
    logical             ::    ifcaneddy   !logical canopy eddy Kz option (default = .FALSE.)
    logical             ::    ifcanphot   !logical canopy photolsyis atten option (default = .FALSE.)
    logical             ::    ifcanbio    !logical canopy biogenic emissions option (default = .FALSE.)
    integer             ::    pai_opt     !integer for PAI values used or calculated (default = 0)
    real(rk)            ::    pai_set     !real value for PAI set values used (default = 4.0)
    integer             ::    lu_opt      !integer for LU type from model mapped to Massman et al. (default = 0/VIIRS)
    integer             ::    z0_opt      !integer for setting first estimate of z0 (default = 0 for Z0_MOD)
    integer             ::    flameh_opt  !Integer for flameh values used or calculated (default = 0)
    real(rk)            ::    flameh_set  !User Set Flame Height (m)
    real(rk)            ::    frp_fac     !FRP tuning factor for flame height calculation (default = 1.0)
    integer             ::    dx_opt      !Integer for dx resolution values used or calculated (default = 0)
    real(rk)            ::    dx_set      !User Set Grid Cell Resolution (m)
    real(rk)            ::    lai_thresh  !User set grid cell LAI threshold to apply canopy conditions (m2/m2)
    real(rk)            ::    frt_thresh  !User set grid cell forest fraction threshold to apply canopy conditions ()
    real(rk)            ::    fch_thresh  !User set grid cell canopy height threshold to apply canopy conditions (m)
    integer             ::    rsl_opt     !RSL option used in model from Rosenzweig et al. 2021 (default = 0, off)
    real(rk)            ::    z0ghc       !ratio of ground roughness length to canopy top height
    real(rk)            ::    lambdars    !Value representing influence of roughness sublayer (nondimensional)
    real(rk)            ::    bio_cce     !MEGAN biogenic emission canopy environment coefficient.
    integer             ::    bio_vert    !MEGAN vertical integration of emissions option (default = 0, off)
    integer             ::    crop_opt    !Set default integer for crop vegtype option from GEDI or user (default = 0)
    real(rk)            ::    crop_set    !Set default value for crop vegtype heights used in model (m) (Default = 3 m)

END MODULE canopy_canopts_mod
