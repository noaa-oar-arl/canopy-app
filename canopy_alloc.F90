
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
! Allocate canopy input variables
!-------------------------------------------------------------------------------

    if(.not.allocated(variables))     allocate(variables(nlat*nlon))
    if(.not.allocated(variables_2d))  allocate(variables_2d(nlon,nlat))

!-------------------------------------------------------------------------------
! Allocate arrays for Internal Canopy Distribution Variables
!-------------------------------------------------------------------------------

    if(.not.allocated(zhc))        allocate(zhc(modlays))
    if(.not.allocated(fafraczInt)) allocate(fafraczInt(modlays))

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Wind Outputs
!-------------------------------------------------------------------------------

    if (ifcanwind .or. ifcanwaf) then
        write(*,*)  'Canopy wind and/or WAF option selected'
        write(*,*)  '-------------------------------'
        if(.not.allocated(canBOT))        allocate(canBOT(modlays))
        if(.not.allocated(canTOP))        allocate(canTOP(modlays))
        if(.not.allocated(canWIND))       allocate(canWIND(nlat*nlon,modlays))
        if(.not.allocated(canWIND_3d))    allocate(canWIND_3d(nlon,nlat,modlays))
        if(.not.allocated(dx))            allocate(dx(nlat*nlon))
        if(.not.allocated(dx_2d))         allocate(dx_2d(nlon,nlat))
        if(.not.allocated(waf))           allocate(waf(nlat*nlon))
        if(.not.allocated(waf_2d))        allocate(waf_2d(nlon,nlat))

    end if

!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Diffusivity Profile Outputs
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        write(*,*)  'Canopy eddy Kz option selected'
        write(*,*)  '-------------------------------'
        if(.not.allocated(Kz))            allocate(Kz(nlat*nlon,modlays))
        if(.not.allocated(Kz_3d))         allocate(Kz_3d(nlon,nlat,modlays))
    end if


!-------------------------------------------------------------------------------
! Allocate arrays for Canopy Photolysis Correction Factor Outputs
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        write(*,*)  'Canopy photolysis option selected'
        write(*,*)  '-------------------------------'
        if(.not.allocated(rjcf))            allocate(rjcf(nlat*nlon,modlays))
        if(.not.allocated(rjcf_3d))         allocate(rjcf_3d(nlon,nlat,modlays))
    end if


END SUBROUTINE canopy_alloc
