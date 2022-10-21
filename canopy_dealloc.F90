
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
    if(allocated(variables))  deallocate(variables)


!-------------------------------------------------------------------------------
! Dellocate arrays for Canopy Distribution
!-------------------------------------------------------------------------------

    if(allocated(zk))         deallocate(zk) !allocated in canopy_readnml
    if(allocated(zhc))        deallocate(zhc)
    if(allocated(fainc))      deallocate(fainc) !allocated in canopy_driver
    if(allocated(fafracz))    deallocate(fafracz) !allocated in canopy_driver
    if(allocated(fafraczInt)) deallocate(fafraczInt)

!-------------------------------------------------------------------------------
! Deallocate arrays for Canopy Wind
!-------------------------------------------------------------------------------

    if (ifcanwind .or. ifcanwaf) then
        if(allocated(canBOT))     deallocate(canBOT)
        if(allocated(canTOP))     deallocate(canTOP)
        if(allocated(canWIND))    deallocate(canWIND)
        if(allocated(dx))         deallocate(dx)
        if(allocated(waf))        deallocate(waf)
    end if

!-------------------------------------------------------------------------------
! Deallocate arrays for Canopy Diffusivity Profile
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        if(allocated(Kz))         deallocate(Kz)
    end if


!-------------------------------------------------------------------------------
! Deallocate arrays for Canopy Photolysis Correction Factor
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        if(allocated(rjcf))         deallocate(rjcf)
    end if

END SUBROUTINE canopy_dealloc
