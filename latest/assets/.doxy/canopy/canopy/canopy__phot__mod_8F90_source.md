

# File canopy\_phot\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_phot\_mod.F90**](canopy__phot__mod_8F90.md)

[Go to the documentation of this file](canopy__phot__mod_8F90.md)


```Fortran



module canopy_phot_mod

    implicit none

contains

    SUBROUTINE canopy_phot( FCLAI, LAI, CLU, COSZEN, RJCF )
        use canopy_const_mod, ONLY: rk     !> constants for canopy models

        REAL(rk),    INTENT( IN )  :: fclai(:)
        REAL(rk),    INTENT( IN )  :: lai
        REAL(rk),    INTENT( IN )  :: clu
        REAL(rk),    INTENT( IN )  :: coszen

        REAL(rk),    INTENT( OUT ) :: rjcf(:)

        rjcf = max(1.0e-10_rk, exp(-1.0_rk*(0.5_rk*(lai*(1.0_rk-fclai))*clu)/max(0.05_rk, coszen)))

    END SUBROUTINE canopy_phot


end module canopy_phot_mod
```


