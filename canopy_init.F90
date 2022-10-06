
SUBROUTINE canopy_init

!-------------------------------------------------------------------------------
! Name:     Initialize Arrays for Canopy Inputs and Outputs
! Purpose:  Initialize arrays for Canopy Inputs and Outputs
! Revised:  03 Oct 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

    USE canopy_canopts_mod
    USE canopy_coord_mod
    USE canopy_canmet_mod
    USE canopy_canvars_mod

    IMPLICIT NONE

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Distribution
!-------------------------------------------------------------------------------

    if(allocated(zhc))        zhc(:)        = 0.0
    if(allocated(fafraczInt)) fafraczInt(:) = 0.0

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Wind
!-------------------------------------------------------------------------------

    if (ifcanwind) then
        if(allocated(canBOT))     canBOT(:)            = 0.0
        if(allocated(canTOP))     canTOP(:)            = 0.0
        if(allocated(canWIND))    canWIND(:,:)         = 0.0
        if(allocated(dx))         dx(:)                = 0.0
        if(allocated(waf))        waf(:)               = 0.0
    end if

!-------------------------------------------------------------------------------
! Initialize arrays for Canopy Diffusivity Profile
!-------------------------------------------------------------------------------

    if (ifcaneddy) then
        if(allocated(Kz))         Kz(:,:) = 0.0
    end if


!-------------------------------------------------------------------------------
! Allocate 3D arrays for Canopy Photolysis Correction Factor
!-------------------------------------------------------------------------------

    if (ifcanphot) then
        if(allocated(rjcf))         rjcf(:,:) = 0.0
    end if


END SUBROUTINE canopy_init
