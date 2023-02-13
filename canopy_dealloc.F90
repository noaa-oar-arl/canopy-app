
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
    if(allocated(variables_2d))  deallocate(variables_2d)

!-------------------------------------------------------------------------------
! Dellocate arrays for Canopy Distribution
!-------------------------------------------------------------------------------

    if(allocated(zk))         deallocate(zk) !allocated in canopy_readnml
    if(allocated(zhc))        deallocate(zhc)
    if(allocated(fainc))      deallocate(fainc) !allocated in canopy_profile
    if(allocated(fafracz))    deallocate(fafracz) !allocated in canopy_profile
    if(allocated(fafraczInt)) deallocate(fafraczInt)

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

    end if

!-------------------------------------------------------------------------------
!  Deallocate NetCDF data structures if used
!-------------------------------------------------------------------------------

    if(allocated(fld2dxy))        deallocate(fld2dxy)
    if(allocated(fld2dxyt))       deallocate(fld2dxyt)
    if(allocated(fld3dxyzt))      deallocate(fld3dxyzt)

END SUBROUTINE canopy_dealloc
