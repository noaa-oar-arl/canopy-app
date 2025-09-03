

# File canopy\_const\_mod.F90

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**canopy\_const\_mod.F90**](canopy__const__mod_8F90.md)

[Go to the documentation of this file](canopy__const__mod_8F90.md)


```Fortran



MODULE canopy_const_mod

    IMPLICIT NONE


    INTEGER, PARAMETER :: rk = selected_real_kind(15, 307)
    REAL(rk),    PARAMETER :: fillreal  = -9.0e20



    REAL(rk),       PARAMETER     :: pi = 3.14159265358979324_rk

    REAL(rk),          PARAMETER     :: pi180 = pi / 180.0_rk



    REAL(rk),          PARAMETER     :: rearth = 6371008.8_rk

    REAL(rk),          PARAMETER     :: siday = 86164.09_rk

    REAL(rk),          PARAMETER     :: grav = 9.80622_rk

    REAL(rk)                         :: dg2m

    REAL(rk),          PARAMETER     :: solcnst = 1373.0




    REAL(rk),          PARAMETER     :: avo = 6.02214076d23

    REAL(rk),          PARAMETER     :: rgasuniv = 8.31446261815324_rk

    REAL(rk),          PARAMETER     :: stdatmpa = 101325.0

    REAL(rk),          PARAMETER     :: stdtemp = 273.15_rk

    REAL(rk),          PARAMETER     :: stfblz = 5.67037442d-8

    REAL(rk),          PARAMETER     :: molvol = 22.4139695_rk




    REAL(rk),          PARAMETER     :: mwair = 28.9628_rk

    REAL(rk),          PARAMETER     :: rdgas = 1.0e3 * rgasuniv / mwair

    REAL(rk),          PARAMETER     :: mwwat = 18.01528_rk

    REAL(rk),          PARAMETER     :: rwvap = 1.0e3 * rgasuniv / mwwat

    REAL(rk),          PARAMETER     :: cpd = 7.0 * rdgas / 2.0

    REAL(rk),          PARAMETER     :: cvd = 5.0 * rdgas / 2.0

    REAL(rk),          PARAMETER     :: cpwvap = 4.0 * rwvap

    REAL(rk),          PARAMETER     :: cvwvap = 3.0 * rwvap

    REAL(rk),          PARAMETER     :: vp0 = 611.29_rk

    REAL(rk),          PARAMETER     :: lv0 = 2.501e6

    REAL(rk),          PARAMETER     :: dlvdt = 2370.0

    REAL(rk),          PARAMETER     :: lf0 = 3.34e5

    REAL(rk),          PARAMETER     :: vonk = 0.4_rk

    REAL(rk),          PARAMETER     :: beta_n = 0.35_rk



    REAL(rk), PARAMETER  :: tau_days  = 5.0_rk

    REAL(rk), PARAMETER  :: tau_hours = 12.0_rk



END MODULE canopy_const_mod
```


