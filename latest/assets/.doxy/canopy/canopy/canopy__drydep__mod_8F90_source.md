

# File canopy\_drydep\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_drydep\_mod.F90**](canopy__drydep__mod_8F90.md)

[Go to the documentation of this file](canopy__drydep__mod_8F90.md)


```Fortran



module canopy_drydep_mod

    implicit none

contains

    SUBROUTINE canopy_gas_drydep_zhang( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        ZK, FCH, TEMPA, PRESSA, &
        RELHUMA, FSUN, PPFD_SUN, PPFD_SHADE, UBAR, &
        SRAD, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !< constants for canopy models
        use canopy_utils_mod,  ONLY: molecdiff,rs_zhang_gas,effhenryslawcoeff,& !< utility functions
            reactivityparam,rbl,rcl,rml

        INTEGER,     INTENT( IN )       :: chemmechgas_opt
        INTEGER,     INTENT( IN )       :: chemmechgas_tot
        REAL(rk),    INTENT( IN )       :: zk(:)
        REAL(rk),    INTENT( IN )       :: fch
        REAL(rk),    INTENT( IN )       :: fsun(:)
        REAL(rk),    INTENT( IN )       :: ppfd_sun(:)
        REAL(rk),    INTENT( IN )       :: ppfd_shade(:)
        REAL(rk),    INTENT( IN )       :: tempa(:)
        REAL(rk),    INTENT( IN )       :: pressa(:)
        REAL(rk),    INTENT( IN )       :: relhuma(:)
        REAL(rk),    INTENT( IN )       :: ubar(:)
        REAL(rk),    INTENT( IN )       :: srad
        REAL(rk),    INTENT( IN )       :: ra
        INTEGER,     INTENT( IN )       :: dep_ind

        REAL(rk),    INTENT( OUT )      :: dep_out(:)

        REAL(rk) ::  ppfd(size(zk))
        REAL(rk) ::  mdiffl(size(zk))
        REAL(rk) ::  rs(size(zk))
        REAL(rk) ::  rb(size(zk))
        REAL(rk) ::  rc(size(zk))
        REAL(rk) ::  rm(size(zk))
        REAL(rk) ::  hstarl
        REAL(rk) ::  f01
        REAL(rk) :: rnum,rden,rlx,vdlx
        INTEGER i

        ppfd = (ppfd_sun*fsun) + (ppfd_shade*(1.0-fsun))

        hstarl  = effhenryslawcoeff(chemmechgas_opt,chemmechgas_tot,dep_ind)
        f01     = reactivityparam(chemmechgas_opt,chemmechgas_tot,dep_ind)
        do i=1, SIZE(zk)
            if (zk(i) .gt. 0.0 .and. zk(i) .le. fch) then           
                mdiffl(i)  = molecdiff(chemmechgas_opt,chemmechgas_tot,dep_ind,tempa(i),pressa(i))
                rs(i)      = rs_zhang_gas(mdiffl(i),tempa(i),pressa(i),ppfd(i),srad,relhuma(i)) 
                rb(i)      = rbl(mdiffl(i), ubar(i)*100.0_rk) 
                rc(i)      = rcl(hstarl, f01)                              
                rm(i)      = rml(hstarl, f01)                              
                rnum = rc(i) * (rs(i) + rm(i))
                rden = rc(i) + 2.0_rk * (rs(i) + rm(i))
                rlx   = rb(i) + (rnum/rden) + ra                           
                vdlx  = 1.0_rk/rlx
                dep_out(i) = vdlx                                          
            else
                rb(i) = 0.0_rk
                rc(i) = 0.0_rk
                rm(i) = 0.0_rk
                rs(i) = 0.0_rk
                dep_out(i) = 0.0_rk
            endif
        end do

    END SUBROUTINE canopy_gas_drydep_zhang

    SUBROUTINE canopy_gas_drydep_soil( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, &
        TEMPSOIL, PRESSA, UBAR, SOCAT, SOTYP, DSOIL, STHETA, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !< constants for canopy models
        use canopy_utils_mod,  ONLY: molecdiff,soilresist,soilrbg !< utility functions

        INTEGER,     INTENT( IN )       :: chemmechgas_opt
        INTEGER,     INTENT( IN )       :: chemmechgas_tot
        REAL(rk),    INTENT( IN )       :: tempsoil
        REAL(rk),    INTENT( IN )       :: pressa
        REAL(rk),    INTENT( IN )       :: ubar
        INTEGER,     INTENT( IN )       :: socat
        INTEGER,     INTENT( IN )       :: sotyp
        REAL(rk),    INTENT( IN )       :: dsoil
        REAL(rk),    INTENT( IN )       :: stheta
        REAL(rk),    INTENT( IN )       :: ra
        INTEGER,     INTENT( IN )       :: dep_ind

        REAL(rk),    INTENT( OUT )      :: dep_out

        real(rk)                        :: mdiffl
        real(rk)                        :: rsoill
        real(rk)                        :: rbg


        mdiffl = molecdiff(chemmechgas_opt,chemmechgas_tot,dep_ind,tempsoil,pressa)

        rsoill = soilresist(mdiffl,socat,sotyp,dsoil,stheta)

        rbg = soilrbg(ubar*100.0_rk)
        dep_out = 1.0_rk/(rbg+rsoill+ra)

        return
    END SUBROUTINE canopy_gas_drydep_soil

    SUBROUTINE canopy_gas_drydep_snow( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, UBAR, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk                       !< constants for canopy models
        use canopy_utils_mod,  ONLY: reactivityparamhno3, soilrbg !< utility functions

        INTEGER,     INTENT( IN )       :: chemmechgas_opt
        INTEGER,     INTENT( IN )       :: chemmechgas_tot
        REAL(rk),    INTENT( IN )       :: ubar
        REAL(rk),    INTENT( IN )       :: ra
        INTEGER,     INTENT( IN )       :: dep_ind

        REAL(rk),    INTENT( OUT )      :: dep_out

        real(rk), parameter             :: ar_0   = 8.0
        real(rk)                        :: ar_l
        real(rk), parameter             :: rsnow0 = 100.0
        real(rk)                        :: rsnowl
        real(rk)                        :: rbg

        ar_l = reactivityparamhno3(chemmechgas_opt,chemmechgas_tot,dep_ind)

        rsnowl = rsnow0 * (ar_0/ar_l)

        rbg = soilrbg(ubar*100.0_rk)

        dep_out = 1.0_rk/(rbg+rsnowl+ra)

        return
    END SUBROUTINE canopy_gas_drydep_snow

    SUBROUTINE canopy_gas_drydep_urban( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, UBAR, TEMP, GAMMA_BUILD, &
        RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk, pi, rgasuniv                       !< constants for canopy models
        use canopy_utils_mod,  ONLY: molarmassgas, soilrbg !< utility functions

        INTEGER,     INTENT( IN )       :: chemmechgas_opt
        INTEGER,     INTENT( IN )       :: chemmechgas_tot
        REAL(rk),    INTENT( IN )       :: ubar
        REAL(rk),    INTENT( IN )       :: temp
        REAL(rk),    INTENT( IN )       :: gamma_build
        REAL(rk),    INTENT( IN )       :: ra
        INTEGER,     INTENT( IN )       :: dep_ind

        REAL(rk),    INTENT( OUT )      :: dep_out

        real(rk)                        :: mmg_l
        real(rk)                        :: cave_l
        real(rk)                        :: rurbanl
        real(rk)                        :: rbg

        mmg_l = molarmassgas(chemmechgas_opt,chemmechgas_tot,dep_ind)

        cave_l = sqrt((8.0_rk*rgasuniv*temp)/(pi*mmg_l))*100.0_rk

        rurbanl = 4.0_rk/(gamma_build*cave_l) 

        rbg = soilrbg(ubar*100.0_rk)

        dep_out = 1.0_rk/(rbg+rurbanl+ra)

        return
    END SUBROUTINE canopy_gas_drydep_urban

    SUBROUTINE canopy_gas_drydep_water( CHEMMECHGAS_OPT,CHEMMECHGAS_TOT, TEMP2, QV2, &
        USTAR, RA, DEP_IND, DEP_OUT)

        use canopy_const_mod,  ONLY: rk, cpd, lv0, dlvdt, stdtemp      !< constants for canopy models
        use canopy_utils_mod,  ONLY: effhenryslawcoeff, lebasmvgas, waterrbw !< utility functions

        INTEGER,     INTENT( IN )       :: chemmechgas_opt
        INTEGER,     INTENT( IN )       :: chemmechgas_tot
        REAL(rk),    INTENT( IN )       :: temp2
        REAL(rk),    INTENT( IN )       :: qv2
        REAL(rk),    INTENT( IN )       :: ustar
        REAL(rk),    INTENT( IN )       :: ra
        INTEGER,     INTENT( IN )       :: dep_ind

        REAL(rk),    INTENT( OUT )      :: dep_out

        real(rk)                        :: hstarl
        real(rk)                        :: ctemp2
        real(rk)                        :: lv
        real(rk)                        :: cp_air
        real(rk)                        :: tw
        real(rk)                        :: lebas_l
        real(rk)                        :: dw
        real(rk)                        :: dw25
        real(rk)                        :: kvisw
        real(rk)                        :: scw_pr_23
        real(rk)                        :: rbw
        real(rk)                        :: rwaterl

        real(rk), Parameter             :: pr         = 0.709
        real(rk), Parameter             :: rt25ink    = 1.0_rk/(stdtemp + 25.0_rk) 
        real(rk), Parameter             :: twothirds  = 2.0_rk / 3.0_rk 
        real(rk), Parameter             :: d3         = 1.38564e-2 

        hstarl  = effhenryslawcoeff(chemmechgas_opt,chemmechgas_tot,dep_ind)

        ctemp2 = temp2 - stdtemp                                      
        lv     = lv0 - dlvdt * ctemp2                                 
        cp_air = cpd * ( 1.0_rk + 0.84_rk * qv2 )                     
        tw     = ( ( 4.71e4 * cp_air / lv ) - 0.870_rk ) + stdtemp    

        hstarl  = hstarl * 0.08205_rk * tw

        lebas_l = lebasmvgas(chemmechgas_opt,chemmechgas_tot,dep_ind)

        dw25 = 13.26e-5 / ( 0.8904_rk**1.14_rk * lebas_l**0.589_rk )
        kvisw = 0.017_rk * exp( -0.025_rk * ( tw - stdtemp ) )
        dw    = dw25 * ( tw * rt25ink ) * ( 0.009025_rk / kvisw )
        scw_pr_23 = ( ( kvisw / dw ) / pr ) ** twothirds

        rwaterl = scw_pr_23 / ( hstarl * d3 * ustar*100.0_rk )

        rbw = waterrbw(d3*ustar*100.0_rk)

        dep_out = 1.0_rk/(rbw+rwaterl+ra)

        return
    END SUBROUTINE canopy_gas_drydep_water


end module canopy_drydep_mod
```


