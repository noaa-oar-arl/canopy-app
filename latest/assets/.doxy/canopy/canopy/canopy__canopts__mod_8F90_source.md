

# File canopy\_canopts\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_canopts\_mod.F90**](canopy__canopts__mod_8F90.md)

[Go to the documentation of this file](canopy__canopts__mod_8F90.md)


```Fortran



MODULE canopy_canopts_mod

!-------------------------------------------------------------------------------
! Name:     Canopy Option Variable Descriptions
! Purpose:  Contains canopy option variable descriptions.
!           03 Oct 2022  Initial Version. (P. C. Campbell)
!-------------------------------------------------------------------------------
    use canopy_const_mod, ONLY: rk
    IMPLICIT NONE

!! .... defines canopy options (read from user namelist)

    integer             ::    infmt_opt

    integer             ::    var3d_opt

    integer             ::    var3d_set

    integer             ::    pavd_opt

    real(rk)            ::    pavd_set

    integer             ::    href_opt

    real(rk)            ::    href_set

    logical             ::    ifcanwind

    logical             ::    ifcanwaf

    logical             ::    ifcaneddy

    logical             ::    ifcanphot

    logical             ::    ifcanbio

    logical             ::    ifcanddepgas

    logical             ::    ifcanaeroddep

    integer             ::    aeroddep_opt

    real(rk)            ::    aeroddep_diam

    real(rk)            ::    aeroddep_rho

    integer             ::    pai_opt

    real(rk)            ::    pai_set

    integer             ::    lu_opt

    integer             ::    z0_opt

    integer             ::    flameh_opt

    integer             ::    flameh_cal

    real(rk)            ::    flameh_set

    real(rk)            ::    frp_fac

    integer             ::    dx_opt

    real(rk)            ::    dx_set

    real(rk)            ::    lai_thresh

    real(rk)            ::    cf_thresh

    real(rk)            ::    ch_thresh

    integer             ::    rsl_opt

    real(rk)            ::    z0ghc

    real(rk)            ::    lambdars
    real(rk)            ::    bio_cce

    integer             ::    biospec_opt

    integer             ::    biovert_opt

    integer             ::    can_opt

    real(rk)            ::    can_chset

    real(rk)            ::    can_cfset

    real(rk)            ::    can_laiset

    integer             ::    ssg_opt

    real(rk)            ::    ssg_chset

    real(rk)            ::    ssg_cfset
    real(rk)            ::    ssg_laiset

    integer             ::    crop_opt

    real(rk)            ::    crop_chset

    real(rk)            ::    crop_cfset

    real(rk)            ::    crop_laiset

    integer             ::    co2_opt

    real(rk)            ::    co2_set

    integer             ::    leafage_opt

    integer             ::    lai_tstep

    integer             ::    loss_opt

    real(rk)            ::    lifetime

    real(rk)            ::    loss_set

    integer             ::    loss_ind

    integer             ::    hist_opt

    integer             ::    soim_opt
    real(rk)            ::    soild1

    real(rk)            ::    soild2

    real(rk)            ::    soild3

    real(rk)            ::    soild4

    integer             ::    aq_opt

    real(rk)            ::    w126_set

    integer             ::    ht_opt

    integer             ::    lt_opt

    integer             ::    hw_opt

    integer             ::    ddepspecgas_opt

    integer             ::    chemmechgas_opt

    integer             ::    chemmechgas_tot

    integer             ::    soilcat_opt

    real(rk)            ::    hyblev1

    real(rk)            ::    snowc_set

    real(rk)            ::    icec_set

    real(rk)            ::    gamma_set

    real(rk)            ::    ramin_set

END MODULE canopy_canopts_mod
```


