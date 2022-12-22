
SUBROUTINE canopy_init

!-------------------------------------------------------------------------------
! Name:     Initialize Arrays for Canopy Inputs and Outputs
! Purpose:  Initialize arrays for Canopy Inputs and Outputs
! Revised:  03 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_const_mod, ONLY: fillreal
    USE canopy_canopts_mod
    USE canopy_coord_mod
    USE canopy_canmet_mod
    USE canopy_canvars_mod

    IMPLICIT NONE

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Distribution
!-------------------------------------------------------------------------------

    if(allocated(zhc))        zhc(:)        = fillreal
    if(allocated(fafraczInt)) fafraczInt(:) = fillreal

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Wind
!-------------------------------------------------------------------------------

    if (ifcanwind .or. ifcanwaf) then
        if(allocated(canBOT))        canBOT(:)            = fillreal
        if(allocated(canTOP))        canTOP(:)            = fillreal
        if(allocated(canWIND))       canWIND(:,:)         = fillreal
        if(allocated(canWIND_3d))    canWIND_3d(:,:,:)    = fillreal
        if(allocated(dx))            dx(:)                = fillreal
        if(allocated(dx_2d))         dx_2d(:,:)           = fillreal
        if(allocated(waf))           waf(:)               = fillreal
        if(allocated(waf_2d))        waf_2d(:,:)          = fillreal
    end if

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Diffusivity Profile
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        if(allocated(Kz))            Kz(:,:)      = fillreal
        if(allocated(Kz_3d))         Kz_3d(:,:,:) = fillreal
    end if


!-------------------------------------------------------------------------------
! Allocate 3D arrays for Canopy Photolysis Correction Factor
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        if(allocated(rjcf))            rjcf(:,:)      = fillreal
        if(allocated(rjcf_3d))         rjcf_3d(:,:,:) = fillreal
    end if


END SUBROUTINE canopy_init
