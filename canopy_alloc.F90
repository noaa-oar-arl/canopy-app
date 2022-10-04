
SUBROUTINE canopy_alloc

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
! Allocate input variables
!-------------------------------------------------------------------------------

    if(.not.allocated(variables))  allocate(variables(nlat*nlon))

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Distribution
!-------------------------------------------------------------------------------

!    if(.not.allocated(zk))         allocate(zk(modlays))
    if(.not.allocated(zhc))        allocate(zhc(modlays))
    if(.not.allocated(fafraczInt)) allocate(fafraczInt(modlays))

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Wind
!-------------------------------------------------------------------------------

    if (ifcanwind) then
        write(*,*)  'Canopy wind/WAF option selected'
        if(.not.allocated(canBOT))     allocate(canBOT(modlays))
        if(.not.allocated(canTOP))     allocate(canTOP(modlays))
        if(.not.allocated(canWIND))    allocate(canWIND(modlays,nlat*nlon))
        if(.not.allocated(dx))         allocate(dx(nlat*nlon))
        if(.not.allocated(waf))        allocate(waf(nlat*nlon))
    end if

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Diffusivity Profile
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        write(*,*)  'Canopy eddy Kz option selected'
        if(.not.allocated(Kz))         allocate(Kz(modlays,nlat*nlon))
    end if


!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Photolysis Correction Factor
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        write(*,*)  'Canopy photolysis option selected'
        if(.not.allocated(rjcf))         allocate(rjcf(modlays,nlat*nlon))
    end if


END SUBROUTINE canopy_alloc
