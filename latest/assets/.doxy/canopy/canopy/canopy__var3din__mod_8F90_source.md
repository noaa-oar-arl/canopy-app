

# File canopy\_var3din\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_var3din\_mod.F90**](canopy__var3din__mod_8F90.md)

[Go to the documentation of this file](canopy__var3din__mod_8F90.md)


```Fortran



module canopy_var3din_mod

    implicit none

contains

    SUBROUTINE canopy_pavd2fafrac ( ZCANMAX_IN, SIGMAU, SIGMA1, FCH, &
        ZHC, PAVD_IN, PAVD_LEVS, FAFRACZINT )
        use canopy_const_mod, ONLY: rk                                           !> constants for canopy models
        use canopy_utils_mod, ONLY: interp_linear1_internal,integratetrapezoid   !> utilities for canopy models

        REAL(rk),    INTENT( IN )  :: fch
        REAL(rk),    INTENT( IN )  :: zcanmax_in
        REAL(rk),    INTENT( IN )  :: sigmau
        REAL(rk),    INTENT( IN )  :: sigma1
        REAL(rk),    INTENT( IN )  :: zhc(:)
        REAL(rk),    INTENT( IN )  :: pavd_in(:)
        REAL(rk),    INTENT( IN )  :: pavd_levs(:)

        REAL(rk),    INTENT( OUT ) :: fafraczint(:)

        INTEGER                     :: i, lev
        REAL(rk)                    :: zk(size(zhc))
        REAL(rk)                    :: pavd_interp(size(zhc))

        REAL(rk)                    :: zcanmax
        REAL(rk), allocatable       :: fainc(:)
        REAL(rk), allocatable       :: fafracz(:)
        REAL(rk)                    :: fatot

        zk = zhc*fch

        pavd_interp = 0.0_rk  !Initialize PAVD_INTERP = 0
        do lev=1, SIZE(pavd_levs) - 1
            do i=2, SIZE(zk)  !loop over only levels ABOVE ground
                if (zk(i) .le.  pavd_levs(1)) then
                    pavd_interp(i)   = pavd_in(1)
                end if
                if (zk(i) .ge.  pavd_levs(lev) .and. zk(i) .le.  pavd_levs(lev+1)) then
                    pavd_interp(i)   = interp_linear1_internal((/ pavd_levs(lev),pavd_levs(lev+1) /), &
                        (/ pavd_in(lev),pavd_in(lev+1) /),zk(i))
                end if
            end do
        end do

        zcanmax             = 0.0_rk

        do i=2, SIZE(zk)
            if (pavd_interp(i) .ge. maxval(pavd_interp) ) then
                zcanmax = zk(i)/fch
                if(zcanmax .gt. 1.0_rk) then !if ZK (at max GEDI PAVD height) > GEDI FCH (inconsistent!)
                    zcanmax = zcanmax_in       !override with Massman Input ZCANMAX
                end if
                exit
            end if
        end do

        if(.not.allocated(fainc))      allocate(fainc(SIZE(zk)))
        if(.not.allocated(fafracz))    allocate(fafracz(SIZE(zk)))

        fainc             = 0.0_rk
        fafracz           = 0.0_rk

        do i=1, SIZE(zk)
            if (zhc(i) >= zcanmax .and. zhc(i) <= 1.0) then
                fainc(i) = exp((-1.0*((zhc(i)-zcanmax)**2.0))/sigmau**2.0)
            else if (zhc(i) >= 0.0 .and. zhc(i) <= zcanmax) then
                fainc(i) = exp((-1.0*((zcanmax-zhc(i))**2.0))/sigma1**2.0)
            end if
        end do
        fatot = integratetrapezoid(zhc,fainc)

        do i=1, SIZE(zk)
            fafracz(i) = fainc(i)/fatot
            fafraczint(i) = integratetrapezoid(zhc(1:i),fafracz(1:i))
        end do

    END SUBROUTINE canopy_pavd2fafrac


end module canopy_var3din_mod
```


