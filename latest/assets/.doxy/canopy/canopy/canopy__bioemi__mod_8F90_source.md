

# File canopy\_bioemi\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_bioemi\_mod.F90**](canopy__bioemi__mod_8F90.md)

[Go to the documentation of this file](canopy__bioemi__mod_8F90.md)


```Fortran



module canopy_bioemi_mod

    implicit none

contains

    SUBROUTINE canopy_bio( ZK, FCLAI, FCH, LAI, FSUN, &
        PPFD_SUN, PPFD_SHADE, TLEAF_SUN, TLEAF_SHADE, PPFD24_SUN, &
        PPFD24_SHADE, TLEAF24_AVE,  &
        PPFD240_SUN, PPFD240_SHADE, &
        TLEAF240_AVE, TKA, DSWRF, TEMP2, LU_OPT, &
        VTYPE, MODRES, CCE, VERT, CO2OPT, CO2SET, &
        LEAFAGEOPT, PASTLAI, CURRENTLAI, TSTEPLAI, &
        LOSSOPT, LOSSSET, LOSSIND, LIFETIME, USTAR, &
        SOIMOPT, SOIM1, SOIM2, SOIM3, SOIM4, SOID1, SOID2, SOID3, &
        SOID4, WILT, AQOPT, W126_SET, W126_REF,  &
        HTOPT, LTOPT, HWOPT, DAILY_MAXT2, DAILY_MINT2, &
        DAILY_MAXWS10, &
        MODLAYS, EMI_IND, EMI_OUT)
        use canopy_const_mod,  ONLY: rk,rgasuniv   !> constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal, &
            get_gamma_co2,get_gamma_leafage, &
            get_gamma_soim, get_gamma_aq, get_gamma_ht, &
            get_gamma_lt, get_gamma_hw, get_canloss_bio
        use canopy_bioparm_mod
        use canopy_tleaf_mod

        REAL(rk),    INTENT( IN )       :: zk(:)
        REAL(rk),    INTENT( IN )       :: fclai(:)
        REAL(rk),    INTENT( IN )       :: fch
        REAL(rk),    INTENT( IN )       :: lai
        REAL(rk),    INTENT( IN )       :: fsun(:)
        REAL(rk),    INTENT( IN )       :: tleaf_sun(:)
        REAL(rk),    INTENT( IN )       :: tleaf_shade(:)
        REAL(rk),    INTENT( IN )       :: ppfd_sun(:)
        REAL(rk),    INTENT( IN )       :: ppfd_shade(:)
        REAL(rk),    INTENT( IN )       :: ppfd24_sun(:)
        REAL(rk),    INTENT( IN )       :: ppfd24_shade(:)
        REAL(rk),    INTENT( IN )       :: tleaf24_ave(:)
        REAL(rk),    INTENT( IN )       :: ppfd240_sun(:)
        REAL(rk),    INTENT( IN )       :: ppfd240_shade(:)
        REAL(rk),    INTENT( IN )       :: tleaf240_ave(:)
        REAL(rk),    INTENT( IN )       :: tka(:)

        REAL(rk),    INTENT( IN )       :: dswrf
        REAL(rk),    INTENT( IN )       :: temp2
        INTEGER,     INTENT( IN )       :: lu_opt
        INTEGER,     INTENT( IN )       :: vtype
        REAL(rk),    INTENT( IN )       :: modres
        REAL(rk),    INTENT( IN )       :: cce
        INTEGER,     INTENT( IN )       :: vert
        INTEGER,     INTENT( IN )       :: co2opt
        REAL(rk),    INTENT( IN )       :: co2set
        INTEGER,     INTENT( IN )       :: soimopt
        REAL(rk),    INTENT( IN )       :: soim1
        REAL(rk),    INTENT( IN )       :: soim2
        REAL(rk),    INTENT( IN )       :: soim3
        REAL(rk),    INTENT( IN )       :: soim4
        REAL(rk),    INTENT( IN )       :: soid1
        REAL(rk),    INTENT( IN )       :: soid2
        REAL(rk),    INTENT( IN )       :: soid3
        REAL(rk),    INTENT( IN )       :: soid4
        REAL(rk),    INTENT( IN )       :: wilt
        INTEGER,     INTENT( IN )       :: aqopt
        REAL(rk),    INTENT( IN )       :: w126_set
        REAL(rk),    INTENT( IN )       :: w126_ref
        INTEGER,     INTENT( IN )       :: htopt
        INTEGER,     INTENT( IN )       :: ltopt
        INTEGER,     INTENT( IN )       :: hwopt
        REAL(rk),    INTENT( IN )       :: daily_maxt2
        REAL(rk),    INTENT( IN )       :: daily_mint2
        REAL(rk),    INTENT( IN )       :: daily_maxws10

        INTEGER,    INTENT( IN )        :: leafageopt
        REAL(rk),    INTENT( IN )       :: pastlai
        REAL(rk),    INTENT( IN )       :: currentlai
        REAL(rk),    INTENT( IN )       :: tsteplai

        INTEGER,     INTENT( IN )       :: modlays
        INTEGER,     INTENT( IN )       :: lossopt
        REAL(rk),    INTENT( IN )       :: lifetime
        REAL(rk),    INTENT( IN )       :: lossset
        INTEGER,     INTENT( IN )       :: lossind

        REAL(rk),    INTENT( IN )       :: ustar
        INTEGER,     INTENT( IN )       :: emi_ind

        REAL(rk),    INTENT( OUT )      :: emi_out(:)


        REAL(rk) :: gammatleaf_sun_ldf_num(size(zk))
        REAL(rk) :: gammatleaf_shade_ldf_num(size(zk))
        REAL(rk) :: gammatleaf_sun_ldf_den(size(zk))
        REAL(rk) :: gammatleaf_shade_ldf_den(size(zk))
        REAL(rk) :: gammatleaf_sun_ldf(size(zk))
        REAL(rk) :: gammatleaf_shade_ldf(size(zk))
        REAL(rk) :: gammatleaf_sun_lif(size(zk))
        REAL(rk) :: gammatleaf_shade_lif(size(zk))
        REAL(rk) :: cp_sun(size(zk))
        REAL(rk) :: cp_shade(size(zk))
        REAL(rk) :: alpha_p_sun(size(zk))
        REAL(rk) :: alpha_p_shade(size(zk))
        REAL(rk) :: gammappfd_sun_ldf(size(zk))
        REAL(rk) :: gammappfd_shade_ldf(size(zk))
        REAL(rk) :: gammatleaf_ppfd_ldf(size(zk))
        REAL(rk) :: gammatleaf_ppfd_lif(size(zk))
        REAL(rk) :: gammatleaf_ppfd_ave(size(zk))
        REAL(rk) :: e_opt(size(zk))
        REAL(rk) :: tleaf_opt(size(zk))
        REAL(rk) :: flai(size(zk))
        REAL(rk) :: vpgwt(size(zk))
        REAL(rk) :: gauss(size(zk))
        REAL(rk) :: ldf
        REAL(rk) :: beta
        REAL(rk) :: ct1
        REAL(rk) :: ceo
        REAL(rk) :: ef
        REAL(rk) :: tabovecanopy

        REAL(rk) :: anew, agro, amat, aold

        REAL(rk) :: caq, taq, dtaq

        REAL(rk) :: cht, tht, dtht

        REAL(rk) :: clt, tlt, dtlt

        REAL(rk) :: chw, thw, dthw

        REAL(rk) :: roota, rootb
        REAL(rk) :: gammasoim

        REAL(rk) :: gammaaq
        REAL(rk) :: gammaht
        REAL(rk) :: gammalt
        REAL(rk) :: gammahw

        REAL(rk) :: gammaco2

        REAL(rk) :: gammaleafage

        REAL(rk) :: canloss_fac

        integer i, layers



        REAL(rk),          PARAMETER     :: ppfd0_sun       =  200.0
        REAL(rk),          PARAMETER     :: ppfd0_shade     =  50.0
        REAL(rk),          PARAMETER     :: ct2             =  230.0_rk


! Calculate maximum normalized emission capacity (E_OPT) and Tleaf at E_OPT
        tleaf_opt = 313.0_rk + (0.6_rk * (tleaf240_ave-297.0_rk)) !Guenther et al. (2012)

! Calculate emission species/plant-dependent mapped emission factors and other important coefficients for gamma terms
        call canopy_biop(emi_ind, lu_opt, vtype, ef, ldf, beta, ct1, ceo, anew, agro, amat, aold, roota, rootb, &
            caq, taq, dtaq, cht, tht, dtht, clt, tlt, dtlt, chw, thw, dthw)

!        print*,'CHT=',CHT,'THT=',THT,'DTHT=',DTHT
!        print*,'CLT=',CLT,'TLT=',TLT,'DTLT=',DTLT
!        print*,'CHW=',CHW,'THW=',THW,'DTHW=',DTHW
        e_opt = ceo * exp(0.05_rk * (tleaf24_ave-297.0_rk)) * exp(0.05_rk * (tleaf240_ave-297.0_rk))

! Calculate gamma (activity) values for average Tleaf (Clifton et al., 2022; based on Guenther et al. 2012)
        gammatleaf_sun_ldf_num = ct2*exp((ct1/(rgasuniv/1000.0))*((1.0/tleaf_opt)-(1.0/tleaf_sun)))
        gammatleaf_sun_ldf_den = (ct2-ct1)*(1.0-exp((ct2/(rgasuniv/1000.0))*((1.0/tleaf_opt)-(1.0/tleaf_sun))))
        gammatleaf_sun_ldf = e_opt*(gammatleaf_sun_ldf_num/gammatleaf_sun_ldf_den)
        gammatleaf_sun_lif = exp(beta*(tka-303.0_rk))

        gammatleaf_shade_ldf_num = ct2*exp((ct1/(rgasuniv/1000.0))*((1.0/tleaf_opt)-(1.0/tleaf_shade)))
        gammatleaf_shade_ldf_den = (ct2-ct1)*(1.0-exp((ct2/(rgasuniv/1000.0))*((1.0/tleaf_opt)-(1.0/tleaf_shade))))
        gammatleaf_shade_ldf = e_opt*(gammatleaf_shade_ldf_num/gammatleaf_shade_ldf_den)
        gammatleaf_shade_lif = exp(beta*(tka-303.0_rk))

! Calculate gamma (activity) values for average PPFD (Clifton et al., 2022; based on Guenther et al. 2012)
        if (dswrf .gt. 0.0_rk) then
            alpha_p_sun = 0.004 - 0.0005*log(ppfd240_sun)
            alpha_p_shade = 0.004 - 0.0005*log(ppfd240_shade)
            cp_sun = 0.0468*(ppfd240_sun**(0.6))*exp(0.005*(ppfd24_sun-ppfd0_sun))
            cp_shade = 0.0468*(ppfd240_shade**(0.6))*exp(0.005*(ppfd24_shade-ppfd0_shade))
            gammappfd_sun_ldf = cp_sun*((alpha_p_sun*ppfd_sun)/sqrt(1.0 + (alpha_p_sun**2.0) * (ppfd_sun**2.0)))
            gammappfd_shade_ldf = cp_shade*((alpha_p_shade*ppfd_shade)/sqrt(1.0 + (alpha_p_shade**2.0) * (ppfd_shade**2.0)))
        else
            gammappfd_sun_ldf = 0.0_rk
            gammappfd_shade_ldf = 0.0_rk
        end if

        gammappfd_sun_ldf = max( gammappfd_sun_ldf, 0.0_rk )
        gammappfd_shade_ldf = max( gammappfd_shade_ldf, 0.0_rk )

        gammatleaf_ppfd_ldf = (gammappfd_sun_ldf*gammatleaf_sun_ldf*fsun) + (gammappfd_shade_ldf*gammatleaf_shade_ldf*(1.0-fsun))
        gammatleaf_ppfd_lif = (gammatleaf_sun_lif*fsun) + (gammatleaf_shade_lif*(1.0-fsun))
        gammatleaf_ppfd_ave = (gammatleaf_ppfd_ldf*ldf) + (gammatleaf_ppfd_lif*(1.0-ldf))
        gammatleaf_ppfd_ave = max( gammatleaf_ppfd_ave, 0.0_rk )

! Get CO2 inhibition factor for isoprene only

        if (emi_ind .eq. 1) then  !Isoprene
            gammaco2 = get_gamma_co2(co2opt,co2set)
        else
            gammaco2 = 1.0_rk
        end if

! Get Soil Moisture Factor
        gammasoim =  get_gamma_soim(soimopt,soim1,soim2,soim3,soim4,soid1,soid2,soid3,soid4,wilt, &
            roota,rootb)

! Get LEAF AGE factor

        tabovecanopy  = temp2   !TEMP2 (above air temp) for TABOVECANOPY
        !do i=1, SIZE(ZK)
        gammaleafage = get_gamma_leafage(leafageopt, pastlai, currentlai, tsteplai, tabovecanopy, anew, agro, amat, aold)
        !end do

! Get AQ Stress Factor
        gammaaq = get_gamma_aq(aqopt, w126_ref, w126_set, caq, taq, dtaq)

! Get HT Stress Factor
        gammaht = get_gamma_ht(htopt, daily_maxt2, cht, tht, dtht)

! Get LT Stress Factor
        gammalt = get_gamma_lt(ltopt, daily_mint2, clt, tlt, dtlt)

! Get HW Stress Factor
        gammahw = get_gamma_hw(hwopt, daily_maxws10, chw, thw, dthw)

! Get canopy loss factor (only used in vertical summing options and empirical formulation and parameters based on isoprene)
! Note:  Allowed for other BVOCs but use caution when applying to compare with above canopy flux observations

        canloss_fac = 1.0_rk  !Initialize
        ! All species
        if (lossind .eq. 0) then
            if (lossopt .eq. 2) then !User set value from NL
                canloss_fac = lossset
            else                !Try and calculate if turned on
                canloss_fac = get_canloss_bio(lossopt,lifetime,ustar,fch)
            end if
        end if

        !Only for a specific biogenic species/indice
        if (lossind .eq. emi_ind) then
            if (lossopt .eq. 2) then !User set value from NL
                canloss_fac = lossset
            else                !Try and calculate if turned on
                canloss_fac = get_canloss_bio(lossopt,lifetime,ustar,fch)
            end if
        end if

! Calculate emissions profile in the canopy
        emi_out = 0.0_rk  ! set initial emissions profile to zero
        flai = 0.0_rk  ! set initial fractional FLAI (LAD) profile to zero

        if (vert .eq. 0) then         !Full 3D leaf-level biogenic emissions (no averaging, summing, or integration)
            do i=1, SIZE(zk)
                if (zk(i) .gt. 0.0 .and. zk(i) .le. fch) then           ! above ground level and at/below canopy top
                    if (i .lt. modlays)  then
                        flai(i) = ((fclai(i+1) - fclai(i)) * lai)/modres    !fractional LAI in each layer converted to LAD (m2 m-3)
                    else
                        flai(i) = flai(modlays-1)
                    end if
                    emi_out(i) = flai(i) * ef * gammatleaf_ppfd_ave(i) * gammaco2 * cce * &
                        gammaleafage * gammasoim * gammaaq * gammaht * gammalt * gammahw ! (ug m-3 hr-1)
                    emi_out(i) = emi_out(i) * 2.7777777777778e-13_rk    !convert emissions output to (kg m-3 s-1)
                end if
            end do
        else if (vert .eq. 1) then       !"MEGANv3-like": Use weighting factors normalized to plant distribution shape (FCLAI)
            !across canopy layers
            layers = min((floor(fch/modres) + 1),modlays)
            do i=1,  SIZE(zk)
                if (zk(i) .gt. 0.0 .and. zk(i) .le. fch) then
                    if (i .lt. modlays)  then
                        flai(i) = ((fclai(i+1) - fclai(i)) * lai)/modres    !fractional LAI in each layer converted to LAD (m2 m-3)
                    else
                        flai(i) = flai(modlays-1)
                    end if
                end if
            end do
            do i=1,  SIZE(zk)
                if (zk(i) .gt. 0.0 .and. zk(i) .le. fch) then
                    vpgwt(i) = (flai(i))/sum(flai(1:layers))
                else
                    vpgwt(i) = 0.0_rk !above canopy
                end if
            end do
            emi_out(SIZE(zk)) = lai * ef * sum(gammatleaf_ppfd_ave(1:layers) * &
                vpgwt(1:layers)) * gammaco2 * cce * gammaleafage * gammasoim * &
                gammaaq * gammaht * gammalt * gammahw * canloss_fac !put into top model layer (ug m-2 hr-1)
            emi_out = emi_out * 2.7777777777778e-13_rk    !convert emissions output to    (kg m-2 s-1)
        else if (vert .eq. 2) then       !"MEGANv3-like": Add weighted sum of activity coefficients using normal distribution
            !across canopy layers using 5 layer numbers directly from MEGANv3
            !--warning: weights are not consistent with FCLAI distribution
            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
            layers = min((floor(fch/modres) + 1),modlays)
            do i=1, SIZE(zk)
                if (zk(i) .gt. fch) then
                    gauss(i) = 0.0
                else if (zk(i) .le. fch .and. zk(i) .gt. fch*(4.0_rk/5.0_rk)) then  !Level 1 - 2
                    gauss(i)   = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                        (/ 0.118464_rk,0.0_rk /),zk(i))
                else if (zk(i) .le. fch*(4.0_rk/5.0_rk) .and. zk(i) .gt. fch*(3.0_rk/5.0_rk)) then  !Level 2 - 3
                    gauss(i)   = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                        (/ 0.239314_rk,0.118464_rk /),zk(i))
                else if (zk(i) .le. fch*(3.0_rk/5.0_rk) .and. zk(i) .gt. fch*(2.0_rk/5.0_rk)) then  !Level 3 - 4
                    gauss(i)   = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                        (/ 0.284444_rk,0.239314_rk /),zk(i))
                else if (zk(i) .le. fch*(2.0_rk/5.0_rk) .and. zk(i) .gt. fch*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                    gauss(i)   = interp_linear1_internal((/ fch*(1.0_rk/5.0_rk),fch*(2.0_rk/5.0_rk) /), &
                        (/ 0.239314_rk,0.284444_rk /),zk(i))
                else if (zk(i) .le. fch*(1.0_rk/5.0_rk) ) then  !Level 4 - Bottom
                    gauss(i)   = interp_linear1_internal((/ zk(1),fch*(1.0_rk/5.0_rk) /), &
                        (/ 0.118464_rk,0.239314_rk /),zk(i))
                end if
            end do

            do i=1, SIZE(zk)
                vpgwt(i) = gauss(i)/sum(gauss(1:layers))
            end do
            emi_out(SIZE(zk)) = lai * ef * sum(gammatleaf_ppfd_ave(1:layers) * &
                vpgwt(1:layers)) * gammaco2 * cce * gammaleafage * gammasoim * &
                gammaaq * gammaht * gammalt * gammahw * canloss_fac    !put into top model layer (ug m-2 hr-1)
            emi_out = emi_out * 2.7777777777778e-13_rk    !convert emissions output to    (kg m-2 s-1)
        else if (vert .eq. 3) then       !"MEGANv3-like": Add weighted sum of activity coefficients equally
            !across canopy layers
            !--warning: weights are not consistent with FCLAI distribution
            !used for biomass distribution used for sunlit/shaded in Gamma TLEAF and GammaPPFD.
            layers = min((floor(fch/modres) + 1),modlays)
            do i=1,  SIZE(zk)
                vpgwt(i) = 1.0_rk/layers
            end do
            emi_out(SIZE(zk)) = lai * ef * sum(gammatleaf_ppfd_ave(1:layers) * &
                vpgwt(1:layers)) * gammaco2 * cce * gammaleafage * gammasoim * &
                gammaaq * gammaht * gammalt * gammahw  * canloss_fac    !put into top model layer (ug m-2 hr-1)
            emi_out = emi_out * 2.7777777777778e-13_rk    !convert emissions output to    (kg m-2 s-1)
        else
            write(*,*)  'Wrong BIOVERT_OPT choice of ', vert, ' in namelist...exiting'
            call exit(2)
        end if

    END SUBROUTINE canopy_bio


end module canopy_bioemi_mod
```


