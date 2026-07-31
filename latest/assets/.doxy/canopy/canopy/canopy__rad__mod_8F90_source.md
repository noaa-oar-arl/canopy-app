

# File canopy\_rad\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_rad\_mod.F90**](canopy__rad__mod_8F90.md)

[Go to the documentation of this file](canopy__rad__mod_8F90.md)


```Fortran



module canopy_rad_mod

    implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_fsun_clu( FCLAI, LAI, CLU, COSZEN, FSUN)

!-----------------------------------------------------------------------

! Description:
!     computes linear interpolation method for PPFD sun/shade in canopy.

! Preconditions:
!     in-canopy height, and model LAI, clumping index, and solar zenith angle

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/23 by PCC
!     Jun 2023 P.C. Campbell: Initial standalone PPFD linear subroutine based on
!                             Silva et al. (2020) exponential curve algorithms
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk     !constants for canopy models
        use canopy_phot_mod

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: FCLAI(:)          ! Input Fractional (z) shapes of the
        ! plant surface distribution (nondimensional), i.e., a Fractional Culmulative LAI
        REAL(RK),    INTENT( IN )       :: LAI               ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )       :: CLU               ! Model input Clumping Index
        REAL(RK),    INTENT( IN )       :: COSZEN            ! Model input Cosine Solar Zenith Angle
        REAL(RK),    INTENT( OUT )      :: FSUN(SIZE(FCLAI)) ! Sunlit/Shaded fraction from photolysis correction factor

!Calculate photolyis shading/correction factor through canopy, i.e., the fraction of sunlit leaves downward through canopy
!  `canopy_phot` gives relative direct beam irradiance,
!  which, multiplied by clumping index, gives sunlit fraction (e.g., Bonan 2019, eq. 14.18)

        call canopy_phot(fclai, lai, clu, coszen, fsun)
        fsun = fsun * clu

    END SUBROUTINE canopy_fsun_clu

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SUBROUTINE canopy_ppfd_exp( ZK, FCH, VBDSF, VDDSF, LAI, FSUN, &
        PPFD_SUN, PPFD_SHADE, PPFD_AVE)

!-----------------------------------------------------------------------

! Description:
!     computes linear interpolation method for PPFD sun/shade in canopy.

! Preconditions:
!     in-canopy height, and model LAI, clumping index, and solar zenith angle

! Subroutines and Functions Called:

! Revision History:
!     Prototype 06/23 by PCC
!     Jun 2023 P.C. Campbell: Initial standalone PPFD linear subroutine based on
!                             Silva et al. (2020) exponential curve algorithms
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        use canopy_const_mod, ONLY: rk     !constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal

! Arguments:
!     IN/OUT
        REAL(RK),    INTENT( IN )       :: ZK(:)                          ! Input model heights (m)
        REAL(RK),    INTENT( IN )       :: FCH                            ! Model input canopy height (m)
        REAL(RK),    INTENT( IN )       :: VBDSF                          ! Model input average surface downward visible flux (W/m2)
        REAL(RK),    INTENT( IN )       :: VDDSF                          ! Model input average surface downward visible flux (W/m2)
        REAL(RK),    INTENT( IN )       :: LAI                            ! Model input total Leaf Area Index
        REAL(RK),    INTENT( IN )       :: FSUN(:)                        ! Sunlit/Shaded fraction from photolysis correction factor
        REAL(RK),    INTENT( OUT )      :: PPFD_SUN(SIZE(ZK))             ! PPFD for sunlit leaves (umol phot/m2 s)
        REAL(RK),    INTENT( OUT )      :: PPFD_SHADE(SIZE(ZK))           ! PPFD for shaded leaves (umol phot/m2 s)
        REAL(RK),    INTENT( OUT )      :: PPFD_AVE(SIZE(ZK))             ! Average PPFD for sunlit and shaded leaves (umol phot/m2 s)

!      LOCAL
        REAL(RK),          PARAMETER     :: CTEMP_1_SUN     =  1.083_rk

        REAL(RK),          PARAMETER     :: CTEMP_2_SUN     =  1.096_rk

        REAL(RK),          PARAMETER     :: CTEMP_3_SUN     =  1.104_rk

        REAL(RK),          PARAMETER     :: CTEMP_4_SUN     =  1.098_rk

        REAL(RK),          PARAMETER     :: CTEMP_5_SUN     =  1.090_rk

        REAL(RK),          PARAMETER     :: DTEMP_1_SUN     =  0.002_rk

        REAL(RK),          PARAMETER     :: DTEMP_2_SUN     =  -0.128_rk

        REAL(RK),          PARAMETER     :: DTEMP_3_SUN     =  -0.298_rk

        REAL(RK),          PARAMETER     :: DTEMP_4_SUN     =  -0.445_rk

        REAL(RK),          PARAMETER     :: DTEMP_5_SUN     =  -0.535_rk

        REAL(RK),          PARAMETER     :: CTEMP_1_SHADE   =  0.871_rk

        REAL(RK),          PARAMETER     :: CTEMP_2_SHADE   =  0.890_rk

        REAL(RK),          PARAMETER     :: CTEMP_3_SHADE   =  0.916_rk

        REAL(RK),          PARAMETER     :: CTEMP_4_SHADE   =  0.941_rk

        REAL(RK),          PARAMETER     :: CTEMP_5_SHADE   =  0.956_rk

        REAL(RK),          PARAMETER     :: DTEMP_1_SHADE   =  0.015_rk

        REAL(RK),          PARAMETER     :: DTEMP_2_SHADE   =  -0.141_rk

        REAL(RK),          PARAMETER     :: DTEMP_3_SHADE   =  -0.368_rk
        REAL(RK),          PARAMETER     :: DTEMP_4_SHADE   =  -0.592_rk

        REAL(RK),          PARAMETER     :: DTEMP_5_SHADE   =  -0.743_rk

        REAL(RK),          PARAMETER     :: FRAC_PAR_VBD        =  0.47_rk
        REAL(RK),          PARAMETER     :: FRAC_PAR_VDD        =  0.57_rk

        REAL(RK) :: PAR_EST

        REAL(RK) :: CTEMP_SUN(SIZE(ZK))

        REAL(RK) :: DTEMP_SUN(SIZE(ZK))

        REAL(RK) :: CTEMP_SHADE(SIZE(ZK))

        REAL(RK) :: DTEMP_SHADE(SIZE(ZK))

        integer i

! Use exponential PPFD model based on Silva et al. (2020) to get approx. sun/shade PPFD
! through canopy
!Citation:
!Silva, S. J., Heald, C. L., and Guenther, A. B.: Development of a reduced-complexity plant canopy
!physics surrogate model for use in chemical transport models: a case study with GEOS-Chem v12.3.0,
!Geosci. Model Dev., 13, 2569–2585, https://doi.org/10.5194/gmd-13-2569-2020, 2020.
        do i=1, SIZE(zk)  !calculate linear change in parameters interpolated to Silva et al. 5 layer canopy regions
            if (zk(i) .gt. fch) then ! above canopy, PPFD_leaf = PPFD_toc (toc=top of canopy)
                ctemp_sun(i)   = 0.0
                dtemp_sun(i)   = 0.0
                ctemp_shade(i) = 0.0
                dtemp_shade(i) = 0.0
            else if (zk(i) .le. fch .and. zk(i) .gt. fch*(4.0_rk/5.0_rk)) then  !Level 1 - 2
                ctemp_sun(i)   = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ ctemp_2_sun,ctemp_1_sun /),zk(i))
                dtemp_sun(i)   = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ dtemp_2_sun,dtemp_1_sun /),zk(i))
                ctemp_shade(i) = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ ctemp_2_shade,ctemp_1_shade /),zk(i))
                dtemp_shade(i) = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ dtemp_2_shade,dtemp_1_shade /),zk(i))
            else if (zk(i) .le. fch*(4.0_rk/5.0_rk) .and. zk(i) .gt. fch*(3.0_rk/5.0_rk)) then  !Level 2 - 3
                ctemp_sun(i)   = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ ctemp_3_sun,ctemp_2_sun /),zk(i))
                dtemp_sun(i)   = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ dtemp_3_sun,dtemp_2_sun /),zk(i))
                ctemp_shade(i) = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ ctemp_3_shade,ctemp_2_shade /),zk(i))
                dtemp_shade(i) = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ dtemp_3_shade,dtemp_2_shade /),zk(i))
            else if (zk(i) .le. fch*(3.0_rk/5.0_rk) .and. zk(i) .gt. fch*(2.0_rk/5.0_rk)) then  !Level 3 - 4
                ctemp_sun(i)   = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ ctemp_4_sun,ctemp_3_sun /),zk(i))
                dtemp_sun(i)   = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ dtemp_4_sun,dtemp_3_sun /),zk(i))
                ctemp_shade(i) = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ ctemp_4_shade,ctemp_3_shade /),zk(i))
                dtemp_shade(i) = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ dtemp_4_shade,dtemp_3_shade /),zk(i))
            else if (zk(i) .le. fch*(2.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                ctemp_sun(i)   = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ ctemp_5_sun,ctemp_4_sun /),zk(i))
                dtemp_sun(i)   = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ dtemp_5_sun,dtemp_4_sun /),zk(i))
                ctemp_shade(i) = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ ctemp_5_shade,ctemp_4_shade /),zk(i))
                dtemp_shade(i) = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ dtemp_5_shade,dtemp_4_shade /),zk(i))
            end if
        end do

        par_est      = (frac_par_vbd * vbdsf) + (frac_par_vdd * vddsf)  !Combining direct and diffuse visible beams --> PAR
        ppfd_sun     = par_est * exp(ctemp_sun + dtemp_sun * lai)  !Silva et al. W/m2 --> umol m-2 s-1
        ppfd_shade   = par_est * exp(ctemp_shade + dtemp_shade * lai)
        ppfd_ave = (ppfd_sun*fsun) + (ppfd_shade*(1.0-fsun)) ! average = sum sun and shade weighted by sunlit fraction

    END SUBROUTINE canopy_ppfd_exp

end module canopy_rad_mod
```


