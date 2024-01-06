
SUBROUTINE canopy_dealloc

!-------------------------------------------------------------------------------
! Name:     Allocate Arrays for Canopy Inputs and Outputs
! Purpose:  Allocate arrays for Canopy Inputs and Outputs
! Revised:  03 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_canopts_mod
    USE canopy_coord_mod
    USE canopy_canmet_mod
    USE canopy_canvars_mod

    IMPLICIT NONE

!-------------------------------------------------------------------------------
! Deallocate input variables
!-------------------------------------------------------------------------------

    if(allocated(variables))     deallocate(variables)
    if(allocated(variables_1d))  deallocate(variables_1d)
    if(allocated(variables_2d))  deallocate(variables_2d)
    if(allocated(variables_3d))  deallocate(variables_3d)
    if(allocated(variables_can))     deallocate(variables_can)

!-------------------------------------------------------------------------------
! Dellocate arrays for Canopy Distribution
!-------------------------------------------------------------------------------

    if(allocated(zk))          deallocate(zk) !allocated in canopy_readnml
    if(allocated(zhc))         deallocate(zhc)
    if(allocated(fainc))       deallocate(fainc) !allocated in canopy_profile
    if(allocated(fafracz))     deallocate(fafracz) !allocated in canopy_profile
    if(allocated(fafraczInt))  deallocate(fafraczInt)
    if(allocated(fsun))        deallocate(fsun)
    if(allocated(tleaf_sun))   deallocate(tleaf_sun)
    if(allocated(tleaf_shade)) deallocate(tleaf_shade)
    if(allocated(tleaf_ave))   deallocate(tleaf_ave)
    if(allocated(ppfd_sun))    deallocate(ppfd_sun)
    if(allocated(ppfd_shade))  deallocate(ppfd_shade)
    if(allocated(ppfd_ave))    deallocate(ppfd_ave)

!-------------------------------------------------------------------------------
! Deallocate arrays for Canopy Wind
!-------------------------------------------------------------------------------

    if (ifcanwind .or. ifcanwaf) then
        if(allocated(canBOT))        deallocate(canBOT)
        if(allocated(canTOP))        deallocate(canTOP)
        if(allocated(canWIND))       deallocate(canWIND)
        if(allocated(canWIND_3d))    deallocate(canWIND_3d)
        if(allocated(dx))            deallocate(dx)
        if(allocated(dx_2d))         deallocate(dx_2d)
        if(allocated(waf))           deallocate(waf)
        if(allocated(waf_2d))        deallocate(waf_2d)
        if(allocated(flameh))        deallocate(flameh)
        if(allocated(flameh_2d))     deallocate(flameh_2d)

    end if

!-------------------------------------------------------------------------------
! Deallocate arrays for Canopy Diffusivity Profile
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        if(allocated(Kz))         deallocate(Kz)
        if(allocated(Kz_3d))      deallocate(Kz_3d)
    end if


!-------------------------------------------------------------------------------
! Deallocate arrays for Canopy Photolysis Correction Factor
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        if(allocated(rjcf))         deallocate(rjcf)
        if(allocated(rjcf_3d))      deallocate(rjcf_3d)
    end if

!-------------------------------------------------------------------------------
! Deallocate arrays for Canopy Biogenic Emissions
!-------------------------------------------------------------------------------

    if (ifcanbio) then
        if(allocated(emi_isop))     deallocate(emi_isop)
        if(allocated(emi_isop_3d))  deallocate(emi_isop_3d)
        if(allocated(emi_myrc))     deallocate(emi_myrc)
        if(allocated(emi_myrc_3d))  deallocate(emi_myrc_3d)
        if(allocated(emi_sabi))     deallocate(emi_sabi)
        if(allocated(emi_sabi_3d))  deallocate(emi_sabi_3d)
        if(allocated(emi_limo))     deallocate(emi_limo)
        if(allocated(emi_limo_3d))  deallocate(emi_limo_3d)
        if(allocated(emi_care))     deallocate(emi_care)
        if(allocated(emi_care_3d))  deallocate(emi_care_3d)
        if(allocated(emi_ocim))     deallocate(emi_ocim)
        if(allocated(emi_ocim_3d))  deallocate(emi_ocim_3d)
        if(allocated(emi_bpin))     deallocate(emi_bpin)
        if(allocated(emi_bpin_3d))  deallocate(emi_bpin_3d)
        if(allocated(emi_apin))     deallocate(emi_apin)
        if(allocated(emi_apin_3d))  deallocate(emi_apin_3d)
        if(allocated(emi_mono))     deallocate(emi_mono)
        if(allocated(emi_mono_3d))  deallocate(emi_mono_3d)
        if(allocated(emi_farn))     deallocate(emi_farn)
        if(allocated(emi_farn_3d))  deallocate(emi_farn_3d)
        if(allocated(emi_cary))     deallocate(emi_cary)
        if(allocated(emi_cary_3d))  deallocate(emi_cary_3d)
        if(allocated(emi_sesq))     deallocate(emi_sesq)
        if(allocated(emi_sesq_3d))  deallocate(emi_sesq_3d)
        if(allocated(emi_mbol))     deallocate(emi_mbol)
        if(allocated(emi_mbol_3d))  deallocate(emi_mbol_3d)
        if(allocated(emi_meth))     deallocate(emi_meth)
        if(allocated(emi_meth_3d))  deallocate(emi_meth_3d)
        if(allocated(emi_acet))     deallocate(emi_acet)
        if(allocated(emi_acet_3d))  deallocate(emi_acet_3d)
        if(allocated(emi_co))       deallocate(emi_co)
        if(allocated(emi_co_3d))    deallocate(emi_co_3d)
        if(allocated(emi_bvoc))     deallocate(emi_bvoc)
        if(allocated(emi_bvoc_3d))  deallocate(emi_bvoc_3d)
        if(allocated(emi_svoc))     deallocate(emi_svoc)
        if(allocated(emi_svoc_3d))  deallocate(emi_svoc_3d)
        if(allocated(emi_ovoc))     deallocate(emi_ovoc)
        if(allocated(emi_ovoc_3d))  deallocate(emi_ovoc_3d)
    end if

!-------------------------------------------------------------------------------
!  Deallocate NetCDF data structures if used
!-------------------------------------------------------------------------------

    if(allocated(fld1dz))        deallocate(fld1dz)
    if(allocated(fld2dxy))        deallocate(fld2dxy)
    if(allocated(fld2dxyt))       deallocate(fld2dxyt)
    if(allocated(fld3dxyzt))      deallocate(fld3dxyzt)

END SUBROUTINE canopy_dealloc
