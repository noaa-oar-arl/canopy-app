MODULE canopy_ncf_io_mod

!-------------------------------------------------------------------------------
! Name:     NetCDF IO
! Purpose:  Contains routines to read NetCDF output.
! Revised:  31 Aug 2004  Original version provided by U. Houston.  (S.-B. Kim)
!           21 Jul 2005  Updated to include error-handling on NetCDF functions
!                        and to conform to NetCDF standard routines for
!                        retrieving data (to eliminate type mismatches).
!                        Updated formatting.  (T. Otte)
!           02 Aug 2005  Changed order of variable declarations in some
!                        subroutines to avoid compile failure on some
!                        machines.  (T. Otte)
!           20 Jun 2006  Removed unused variables.  Changed local variables
!                        FILE to FILENAME and DATA to DATAOUT to avoid
!                        conflicts with F90 keywords.  (T. Otte)
!           19 Apr 2007  Added new routine GET_VAR_REAL2_CDF to read "real"
!                        scalars, as needed for WRFv2.2.  Added new routine
!                        GET_DIM_INT_CDF to retrieve netCDF dimensions.
!                        Changed internal error handling so that errors are
!                        passed back using a non-zero RCODE for dimension
!                        mismatches in addition to netCDF errors.  (T. Otte)
!           12 Feb 2010  Removed unused variable ID_TIME from subroutine
!                        GET_GL_ATT_INT_CDF.  (T. Otte)
!           19 Mar 2010  Removed routines GET_DIMS_CDF, GET_DIM_INT_CDF,
!                        GET_GL_ATT_INT_CDF, GET_GL_ATT_REAL_CDF,
!                        GET_GL_ATT_TEXT_CDF, and GET_DIM_ATT_INT_CDF.
!                        Removed file open and close functions from all
!                        remaining routines, and changed input argument
!                        from FILENAME to CDFID.  (T. Otte)
!           31 Aug 2011  Updated netCDF to F90.  Removed GET_TIMES_CDF.
!                        Changed F77 character declarations to F90 standard.
!                        (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           02 Feb 2018  Added new routine GET_VAR_3D_INT_CDF.  (T. Spero)
!           22 Jun 2018  Changed module name from WRF_NETCDF to NETCDF_IO.
!                        (T. Spero)
!           19 Dec 2022  Brought module to the Canopy-App (P.C. Campbell)
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_3d_real_cdf (cdfid, var, dum3d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum3d    ( : , : , : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER                          :: nz
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var


        nx = SIZE(dum3d,1)
        ny = SIZE(dum3d,2)
        nz = SIZE(dum3d,3)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, dum3d, start=(/1,1,1,it/),  &
            count=(/nx,ny,nz,1/))

        IF ( rcode /= nf90_noerr ) then
            print*,'read error ',cdfid,var
            print*,'nx,ny,nz=',nx,ny,nz
        endif
    END SUBROUTINE get_var_3d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_3d_int_cdf (cdfid, var, idum3d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: idum3d    ( : , : , : )
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER                          :: nz
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(idum3d,1)
        ny = SIZE(idum3d,2)
        nz = SIZE(idum3d,3)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, idum3d, start=(/1,1,1,it/),  &
            count=(/nx,ny,nz,1/))

    END SUBROUTINE get_var_3d_int_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_2d_real_cdf (cdfid, var, dum2d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum2d    ( : , : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(dum2d,1)
        ny = SIZE(dum2d,2)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) then
            print*,'can not find variable ',trim(var)
            RETURN
        endif
        rcode = nf90_get_var (cdfid, id_data, dum2d, start=(/1,1,it/),  &
            count=(/nx,ny,1/))
        IF ( rcode /= nf90_noerr ) print*,'read error for ',trim(var),&
            ' id_data,it,nx,ny=',id_data,it,nx,ny

    END SUBROUTINE get_var_2d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_2d_int_cdf (cdfid, var, idum2d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: idum2d   ( : , : )
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER                          :: ny
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(idum2d,1)
        ny = SIZE(idum2d,2)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, idum2d, start=(/1,1,it/),  &
            count=(/nx,ny,1/))

    END SUBROUTINE get_var_2d_int_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_1d_real_cdf (cdfid, var, dum1d, it, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL(rk),          INTENT(OUT)   :: dum1d    ( : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(dum1d)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, dum1d, start=(/1,it/),  &
            count=(/nx,1/))

    END SUBROUTINE get_var_1d_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_1d_int_cdf (cdfid, var, idum1d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: idum1d   ( : )
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        nx = SIZE(idum1d)

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, idum1d, start=(/1,it/),  &
            count=(/nx,1/))

    END SUBROUTINE get_var_1d_int_cdf

!-------
    SUBROUTINE get_var_1d_double_cdf (cdfid, var, dum1d, it, rcode)

        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        REAL,  INTENT(OUT)               :: dum1d    ( : )
        INTEGER                          :: id_data
        INTEGER,           INTENT(IN)    :: it
        INTEGER                          :: nx
        INTEGER,           INTENT(OUT)   :: rcode
        CHARACTER(LEN=*),  INTENT(IN)    :: var
        double precision, allocatable  :: dbtmp(:)

        nx = SIZE(dum1d)
        allocate(dbtmp(nx))

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, dbtmp, start=(/1,it/),  &
            count=(/nx,1/))
        dum1d(:)=sngl(dbtmp(:))
    END SUBROUTINE get_var_1d_double_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

    SUBROUTINE get_var_real_cdf (cdfid, var, scalar, rcode)

        USE canopy_const_mod, ONLY: rk
        USE netcdf

        IMPLICIT NONE

        INTEGER,           INTENT(IN)    :: cdfid
        INTEGER                          :: id_data
        INTEGER,           INTENT(OUT)   :: rcode
        REAL(rk),          INTENT(OUT)   :: scalar
        CHARACTER(LEN=*),  INTENT(IN)    :: var

        rcode = nf90_inq_varid (cdfid, var, id_data)
        IF ( rcode /= nf90_noerr ) RETURN

        rcode = nf90_get_var (cdfid, id_data, scalar)

    END SUBROUTINE get_var_real_cdf

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

END MODULE canopy_ncf_io_mod
