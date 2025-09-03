

# File canopy\_tleaf\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_tleaf\_mod.F90**](canopy__tleaf__mod_8F90.md)

[Go to the documentation of this file](canopy__tleaf__mod_8F90.md)


```Fortran



module canopy_tleaf_mod

    implicit none

contains

    SUBROUTINE canopy_tleaf_lin( ZK, FCH, TEMP2, FSUN, &
        TLEAF_SUN, TLEAF_SHADE, TLEAF_AVE)

        use canopy_const_mod, ONLY: rk     !< constants for canopy models
        use canopy_utils_mod,  ONLY: interp_linear1_internal !< linear interpolation utility

        REAL(rk),    INTENT( IN )       :: zk(:)
        REAL(rk),    INTENT( IN )       :: fch
        REAL(rk),    INTENT( IN )       :: temp2
        REAL(rk),    INTENT( IN )       :: fsun(:)

        REAL(rk),    INTENT( OUT )      :: tleaf_sun(size(zk))
        REAL(rk),    INTENT( OUT )      :: tleaf_shade(size(zk))
        REAL(rk),    INTENT( OUT )      :: tleaf_ave(size(zk))

        REAL(rk),          PARAMETER     :: atemp_1_sun     =  -13.891_rk
        REAL(rk),          PARAMETER     :: atemp_2_sun     =  -12.322_rk
        REAL(rk),          PARAMETER     :: atemp_3_sun     =  -1.032_rk
        REAL(rk),          PARAMETER     :: atemp_4_sun     =  -5.172_rk
        REAL(rk),          PARAMETER     :: atemp_5_sun     =  -5.589_rk
        REAL(rk),          PARAMETER     :: btemp_1_sun     =  1.064_rk
        REAL(rk),          PARAMETER     :: btemp_2_sun     =  1.057_rk
        REAL(rk),          PARAMETER     :: btemp_3_sun     =  1.031_rk
        REAL(rk),          PARAMETER     :: btemp_4_sun     =  1.050_rk
        REAL(rk),          PARAMETER     :: btemp_5_sun     =  1.051_rk
        REAL(rk),          PARAMETER     :: atemp_1_shade   =  -12.846_rk
        REAL(rk),          PARAMETER     :: atemp_2_shade   =  -11.343_rk
        REAL(rk),          PARAMETER     :: atemp_3_shade   =  -1.068_rk
        REAL(rk),          PARAMETER     :: atemp_4_shade   =  -5.551_rk
        REAL(rk),          PARAMETER     :: atemp_5_shade   =  -5.955_rk
        REAL(rk),          PARAMETER     :: btemp_1_shade   =  1.060_rk
        REAL(rk),          PARAMETER     :: btemp_2_shade   =  1.053_rk
        REAL(rk),          PARAMETER     :: btemp_3_shade   =  1.031_rk
        REAL(rk),          PARAMETER     :: btemp_4_shade   =  1.051_rk
        REAL(rk),          PARAMETER     :: btemp_5_shade   =  1.053_rk

        REAL(rk) :: atemp_sun(size(zk))
        REAL(rk) :: btemp_sun(size(zk))
        REAL(rk) :: atemp_shade(size(zk))
        REAL(rk) :: btemp_shade(size(zk))

        integer i

        do i=1, SIZE(zk)  
            if (zk(i) .gt. fch) then 
                atemp_sun(i)   = 0.0
                btemp_sun(i)   = 1.0
                atemp_shade(i) = 0.0
                btemp_shade(i) = 1.0
            else if (zk(i) .le. fch .and. zk(i) .gt. fch*(4.0_rk/5.0_rk)) then  
                atemp_sun(i)   = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ atemp_2_sun,atemp_1_sun /),zk(i))
                btemp_sun(i)   = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ btemp_2_sun,btemp_1_sun /),zk(i))
                atemp_shade(i) = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ atemp_2_shade,atemp_1_shade /),zk(i))
                btemp_shade(i) = interp_linear1_internal((/ fch*(4.0_rk/5.0_rk),fch /), &
                    (/ btemp_2_shade,btemp_1_shade /),zk(i))
            else if (zk(i) .le. fch*(4.0_rk/5.0_rk) .and. zk(i) .gt. fch*(3.0_rk/5.0_rk)) then  
                atemp_sun(i)   = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ atemp_3_sun,atemp_2_sun /),zk(i))
                btemp_sun(i)   = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ btemp_3_sun,btemp_2_sun /),zk(i))
                atemp_shade(i) = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ atemp_3_shade,atemp_2_shade /),zk(i))
                btemp_shade(i) = interp_linear1_internal((/ fch*(3.0_rk/5.0_rk),fch*(4.0_rk/5.0_rk) /), &
                    (/ btemp_3_shade,btemp_2_shade /),zk(i))
            else if (zk(i) .le. fch*(3.0_rk/5.0_rk) .and. zk(i) .gt. fch*(2.0_rk/5.0_rk)) then  
                atemp_sun(i)   = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ atemp_4_sun,atemp_3_sun /),zk(i))
                btemp_sun(i)   = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ btemp_4_sun,btemp_3_sun /),zk(i))
                atemp_shade(i) = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ atemp_4_shade,atemp_3_shade /),zk(i))
                btemp_shade(i) = interp_linear1_internal((/ fch*(2.0_rk/5.0_rk),fch*(3.0_rk/5.0_rk) /), &
                    (/ btemp_4_shade,btemp_3_shade /),zk(i))
            else if (zk(i) .le. fch*(2.0_rk/5.0_rk) ) then  
                atemp_sun(i)   = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ atemp_5_sun,atemp_4_sun /),zk(i))
                btemp_sun(i)   = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ btemp_5_sun,btemp_4_sun /),zk(i))
                atemp_shade(i) = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ atemp_5_shade,atemp_4_shade /),zk(i))
                btemp_shade(i) = interp_linear1_internal((/ zk(1),fch*(2.0_rk/5.0_rk) /), &
                    (/ btemp_5_shade,btemp_4_shade /),zk(i))
            end if
        end do

        tleaf_sun   = atemp_sun + (btemp_sun*temp2)
        tleaf_shade = atemp_shade + (btemp_shade*temp2)
        tleaf_ave = (tleaf_sun*fsun) + (tleaf_shade*(1.0-fsun))

    END SUBROUTINE canopy_tleaf_lin


end module canopy_tleaf_mod
```


