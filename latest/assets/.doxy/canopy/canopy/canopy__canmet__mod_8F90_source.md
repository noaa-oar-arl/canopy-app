

# File canopy\_canmet\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_canmet\_mod.F90**](canopy__canmet__mod_8F90.md)

[Go to the documentation of this file](canopy__canmet__mod_8F90.md)


```Fortran



MODULE canopy_canmet_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Met/Sfc Input Variable Descriptions
! Purpose:  Contains canopy met and sfc input  variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk

    IMPLICIT NONE
!! .... defines canopy options (read from user namelist)


    TYPE :: variable_type
        real(rk)   :: lat

        real(rk)   :: lon

        real(rk)   :: ch

        real(rk)   :: ugrd10m

        real(rk)   :: vgrd10m

        real(rk)   :: clu

        real(rk)   :: lai

        integer    :: vtype

        real(rk)   :: canfrac

        real(rk)   :: fricv

        real(rk)   :: csz

        real(rk)   :: sfcr

        real(rk)   :: mol

        real(rk)   :: frp

        real(rk)   :: href

        integer    :: sotyp

        real(rk)   :: pressfc

        real(rk)   :: dswrf

        real(rk)   :: shtfl

        real(rk)   :: tmpsfc

        real(rk)   :: tmp2m

        real(rk)   :: spfh2m

        real(rk)   :: hpbl

        real(rk)   :: prate_ave

        real(rk)   :: soilw1

        real(rk)   :: soilw2

        real(rk)   :: soilw3

        real(rk)   :: soilw4

        real(rk)   :: wilt

        real(rk)   :: ozone_w126

        real(rk)   :: soilt1

        real(rk)   :: soilt2

        real(rk)   :: soilt3

        real(rk)   :: soilt4

        real(rk)   :: tmp_hyblev1

        real(rk)   :: snowc_ave

        real(rk)   :: icec

        real(rk)   :: vbdsf_ave

        real(rk)   :: vddsf_ave

    end TYPE variable_type

    type(variable_type), allocatable :: variables( : ), variables_2d( : , :)

    TYPE :: variable_type_1d
        real(rk)   :: lev
    end TYPE variable_type_1d

    type(variable_type_1d), allocatable :: variables_1d( : )

    TYPE :: variable_type_3d
        real(rk)   :: pavd
    end TYPE variable_type_3d

    type(variable_type_3d), allocatable :: variables_3d( : , : , :)

    TYPE :: variable_type_can
        real(rk)   :: lat

        real(rk)   :: lon

        real(rk)   :: lev01

        real(rk)   :: pavd01

        real(rk)   :: lev02
        real(rk)   :: pavd02
        real(rk)   :: lev03
        real(rk)   :: pavd03
        real(rk)   :: lev04
        real(rk)   :: pavd04
        real(rk)   :: lev05
        real(rk)   :: pavd05
        real(rk)   :: lev06
        real(rk)   :: pavd06
        real(rk)   :: lev07
        real(rk)   :: pavd07
        real(rk)   :: lev08
        real(rk)   :: pavd08
        real(rk)   :: lev09
        real(rk)   :: pavd09
        real(rk)   :: lev10
        real(rk)   :: pavd10
        real(rk)   :: lev11
        real(rk)   :: pavd11
        real(rk)   :: lev12
        real(rk)   :: pavd12
        real(rk)   :: lev13
        real(rk)   :: pavd13
        real(rk)   :: lev14
        real(rk)   :: pavd14
    end TYPE variable_type_can

    type(variable_type_can), allocatable :: variables_can( : )



    real(rk)       ::    latref

    real(rk)       ::    lonref

    real(rk)       ::    hcmref

    real(rk)       ::    uref

    real(rk)       ::    vref

    real(rk)       ::    ubzref

    real(rk)       ::    cluref

    real(rk)       ::    lairef

    integer        ::    vtyperef

    real(rk)       ::    canfracref

    real(rk)       ::    ustref

    real(rk)       ::    cszref

    real(rk)       ::    z0ref

    real(rk)       ::    molref

    real(rk)       ::    frpref

    real(rk)       ::    hgtref

    integer        ::    sotypref

    real(rk)       ::    pressfcref

    real(rk)       ::    dswrfref

    real(rk)       ::    shtflref

    real(rk)       ::    tmpsfcref

    real(rk)       ::    tmp2mref

    real(rk)       ::    spfh2mref

    real(rk)       ::    hpblref

    real(rk)       ::    prate_averef

    real(rk)       ::    soilw1ref

    real(rk)       ::    soilw2ref

    real(rk)       ::    soilw3ref

    real(rk)       ::    soilw4ref

    real(rk)       ::    wiltref

    real(rk)       ::    ozone_w126ref

    real(rk)       ::    soilt1ref

    real(rk)       ::    soilt2ref

    real(rk)       ::    soilt3ref

    real(rk)       ::    soilt4ref

    real(rk)       ::    tmp_hyblev1ref

    real(rk)       ::    snowc_averef

    real(rk)       ::    icec_averef

    real(rk)       ::    vbdsf_averef

    real(rk)       ::    vddsf_averef

!    real(rk)       ::    lev01ref, lev02ref, lev03ref, lev04ref, lev05ref, & !Input canopy profile levels
!                         lev06ref, lev07ref, lev08ref, lev09ref, lev10ref, &
!                         lev11ref, lev12ref, lev13ref, lev14ref
!    real(rk)       ::    pavd01ref, pavd02ref, pavd03ref, pavd04ref, pavd05ref, & !Input canopy PAVD profile
!                         pavd06ref, pavd07ref, pavd08ref, pavd09ref, pavd10ref, &
!                         pavd11ref, pavd12ref, pavd13ref, pavd14ref
    real(rk), allocatable   :: pavdref ( : ), pavd_arr ( : )

    real(rk), allocatable   :: levref ( : ), lev_arr  ( : )


END MODULE canopy_canmet_mod
```


