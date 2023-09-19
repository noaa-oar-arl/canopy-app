module canopy_utils_mod

    use canopy_const_mod, ONLY: pi, rk, rearth    !constants for canopy models

    implicit none

    private
    public IntegrateTrapezoid,interp_linear1_internal,CalcPAI, &
        CalcDX,CalcFlameH,GET_GAMMA_CO2,GET_GAMMA_LEAFAGE

contains

    function IntegrateTrapezoid(x, y)
        !! Calculates the integral of an array y with respect to x using the trapezoid
        !! approximation. Note that the mesh spacing of x does not have to be uniform.
        real(rk), intent(in)  :: x(:)                !! Variable x
        real(rk), intent(in)  :: y(size(x))          !! Function y(x)
        real(rk)              :: IntegrateTrapezoid  !! Integral ∫y(x)·dx
        ! Integrate using the trapezoidal rule
        associate(n => size(x))
            IntegrateTrapezoid = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
        end associate
    end function
!-------------------------------------------------------------------------------------

    function interp_linear1_internal(x,y,xout) result(yout)
        !! Interpolates for the y value at the desired x value,
        !! given x and y values around the desired point.

        implicit none

        real(rk), intent(IN)  :: x(2), y(2), xout
        real(rk) :: yout
        real(rk) :: alph

        if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
            write(*,*) "interp1: xout < x0 or xout > x1 !"
            write(*,*) "xout = ",xout
            write(*,*) "x0   = ",x(1)
            write(*,*) "x1   = ",x(2)
            stop
        end if

        alph = (xout - x(1)) / (x(2) - x(1))
        yout = y(1) + alph*(y(2) - y(1))

        return

    end function interp_linear1_internal
!--------------------------------------------------------------------------------------

    function CalcPAI(fch, ffrac)
        !! Calculates the Plant Area Index as a function of canopy height and canopy/
        !! forest fraction (Based on Eq. 19 of Massman et al., 2017).

        !!  W.J. Massman, J.M. Forthofer, and M.A. Finney. An improved
        !!  canopy wind model for predicting wind adjustment factors
        !!  and wildland fire behavior. Canadian Journal of Forest Research.
        !!  47(5): 594-603. https://doi.org/10.1139/cjfr-2016-0354

        !! Assume Canopy Cover Fraction, C,  = FFRAC
        !! Assume Canopy Crown Ratio, F, = CC/3.0 = FFRAC/3.0 (Eq. 9 in  Andrews, 2012).
        !! Andrews, P.L. 2012. Modeling wind adjustment factor and midflame wind speed
        !!  for Rothermel’s surface fire spread model. USDA For. Serv. Gen. Tech. Rep. RMRS-GTR-266.


        real(rk), intent(in)  :: fch          !! Input Grid cell canopy height (m
        real(rk), intent(in)  :: ffrac        !! Input Grid cell forest fraction
        real(rk)              :: CalcPAI      !! Calculated Plant area index (PAI)

        CalcPAI=( (fch*(ffrac/3.0_rk)*10.6955_rk) / (2.0_rk * pi) ) * ffrac !Massman PAI calculation (Eq. 19)

    end function
    !--------------------------------------------------------------------------------------

    real(rk) function CalcDX(lat, dlon) result(dx)
        !! Compute the zonal distance, dx, corresponding to longitude increment `dlon`.

        real(rk), intent(in) :: lat   !! Latitude (degrees)
        real(rk), intent(in) :: dlon  !! Longitude increment (degrees)
        real(rk) :: lat_rad, dlon_rad

        lat_rad = lat * pi / 180._rk
        dlon_rad = dlon * pi / 180._rk

        dx = rearth * cos(lat_rad) * dlon_rad

    end function
    !--------------------------------------------------------------------------------------

    ! real(rk) function CalcGCDist(lat1, lat2, lon1, lon2) result(d)
    !     !! Compute great-circle distance between two points using the spherical law of cosines formula.

    !     real(rk), intent(in)  :: lat1,lat2                  !! Two model latitudes
    !     real(rk), intent(in)  :: lon1,lon2                  !! Two model longitudes
    !     real(rk) :: lat_rad1, lat_rad2, lon_rad1, lon_rad2  !! radians

    !     lat_rad1 = lat1/(180.0_rk/pi)
    !     lon_rad1 = lon1/(180.0_rk/pi)
    !     lat_rad2 = lat2/(180.0_rk/pi)
    !     lon_rad2 = lon2/(180.0_rk/pi)

    !     d = rearth*acos( &
    !         sin(lat_rad1)*sin(lat_rad2) &
    !         + cos(lat_rad1)*cos(lat_rad2)*cos(lon_rad2-lon_rad1) &
    !     )

    ! end function
    ! !--------------------------------------------------------------------------------------

    function CalcFlameH(frp, dx)
        !! Approximates the Flame Height as a function of FRP intensity and grid cell distance (dx)
        !! forest fraction (Based on Byram 1959).

        !!  Byram, GM (1959). Combustion of Forest Fuels. In Forest Fire: Control and Use.
        !!  (Ed. KP David) pp. 61-89.  McGraw Hill, New York, NY

        !! Assume Flame Length = Flame Height under calm winds

        real(rk), intent(in)  :: frp          !! Input Grid cell Fire Radiative Power (MW/cell)
        real(rk), intent(in)  :: dx           !! Input Grid cell length (m)
        real(rk)              :: CalcFlameH   !! Calculated Plant area index (PAI)

        CalcFlameH=0.0775_rk*((frp*1000.0_rk)/dx)**0.46_rk  !Byram flameh calculation as function of FRP

    end function

    real(rk) function GET_GAMMA_CO2(co2_opt, co2_set)  result( GAMMA_CO2 )
        ! !IROUTINE: get_gamma_co2
        !
        ! !DESCRIPTION: Function GET\_GAMMA\_CO2 computes the CO2 activity factor
        !  associated with CO2 inhibition of isoprene emission. Called from
        !  GET\_MEGAN\_EMISSIONS only.
        !\\
        !\\
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        integer,  INTENT(IN) :: co2_opt       ! Option for co2 inhibition calculation
        ! 0=Possell & Hewitt (2011);
        ! 1=Wilkinson et al. (2009)
        ! >1=off
        real(rk), INTENT(IN) :: co2_set       ! User set atmospheric CO2 conc [ppmv]
        !
        ! !RETURN VALUE:
!        REAL(rk)             :: GAMMA_CO2  ! CO2 activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        REAL(rk)             :: CO2i       ! Intercellular CO2 conc [ppmv]
        REAL(rk)             :: ISMAXi     ! Asymptote for intercellular CO2
        REAL(rk)             :: HEXPi      ! Exponent for intercellular CO2
        REAL(rk)             :: CSTARi     ! Scaling coef for intercellular CO2
        REAL(rk)             :: ISMAXa     ! Asymptote for atmospheric CO2
        REAL(rk)             :: HEXPa      ! Exponent for atmospheric CO2
        REAL(rk)             :: CSTARa     ! Scaling coef for atmospheric CO2
        !
        ! !REMARKS:
        !  References:
        !  ============================================================================
        !  (1 ) Heald, C. L., Wilkinson, M. J., Monson, R. K., Alo, C. A.,
        !       Wang, G. L., and Guenther, A.: Response of isoprene emission
        !       to ambient co(2) changes and implications for global budgets,
        !       Global Change Biology, 15, 1127-1140, 2009.
        !  (2 ) Wilkinson, M. J., Monson, R. K., Trahan, N., Lee, S., Brown, E.,
        !       Jackson, R. B., Polley, H. W., Fay, P. A., and Fall, R.: Leaf
        !       isoprene emission rate as a function of atmospheric CO2
        !       concentration, Global Change Biology, 15, 1189-1200, 2009.
        !  (3 ) Possell, M., and Hewitt, C. N.: Isoprene emissions from plants
        !       are mediated by atmospheric co2 concentrations, Global Change
        !       Biology, 17, 1595-1610, 2011.

        ! !REVISION HISTORY:
        !  (1 ) Implemented in the standard code by A. Tai (Jun 2012).
        !  See https://github.com/geoschem/hemco for complete history
        !EOP
        !------------------------------------------------------------------------------
        !BOC

        !-----------------------
        ! Compute GAMMA_CO2
        !-----------------------

        !----------------------------------------------------------
        ! Choose between two alternative CO2 inhibition schemes
        !----------------------------------------------------------

        IF ( co2_opt .eq. 0 ) THEN

            ! Use empirical relationship of Possell & Hewitt (2011):
            ! Empirical relationship of Possell & Hewitt (2011) based on nine
            ! experimental studies including Wilkinson et al. (2009).

            GAMMA_CO2 = 8.9406_rk / ( 1.0_rk + 8.9406_rk * 0.0024_rk * co2_set )


        ELSEIF ( co2_opt .eq. 1 ) THEN

            ! Use parameterization of Wilkinson et al. (2009):
            ! Semi-process-based parameterization of Wilkinson et al. (2009),
            ! taking into account of sensitivity to intercellular CO2
            ! fluctuation, which is here set as a constant fraction of
            ! atmospheric CO2. This is especially recommended for sub-ambient
            ! CO2 concentrations::


            ! Parameters for intercellular CO2 using linear interpolation:
            IF ( co2_set <= 600.0_rk ) THEN
                ISMAXi = 1.036_rk  - (1.036_rk - 1.072_rk) / &
                    (600.0_rk - 400.0_rk) * (600.0_rk - co2_set)
                HEXPi  = 2.0125_rk - (2.0125_rk - 1.7000_rk) / &
                    (600.0_rk - 400.0_rk) * (600.0_rk - co2_set)
                CSTARi = 1150.0_rk - (1150.0_rk - 1218.0_rk) / &
                    (600.0_rk - 400.0_rk) * (600.0_rk - co2_set)
            ELSEIF ( co2_set > 600.0_rk .AND. co2_set < 800.0_rk ) THEN
                ISMAXi = 1.046_rk  - (1.046_rk - 1.036_rk) / &
                    (800.0_rk - 600.0_rk) * (800.0_rk - co2_set)
                HEXPi  = 1.5380_rk - (1.5380_rk - 2.0125_rk) / &
                    (800.0_rk - 600.0_rk) * (800.0_rk - co2_set)
                CSTARi = 2025.0_rk - (2025.0_rk - 1150.0_rk) / &
                    (800.0_rk - 600.0_rk) * (800.0_rk - co2_set)
            ELSE
                ISMAXi = 1.014_rk - (1.014_rk - 1.046_rk) / &
                    (1200.0_rk - 800.0_rk) * (1200.0_rk - co2_set)
                HEXPi  = 2.8610_rk - (2.8610_rk - 1.5380_rk) / &
                    (1200.0_rk - 800.0_rk) * (1200.0_rk - co2_set)
                CSTARi = 1525.0_rk - (1525.0_rk - 2025.0_rk) / &
                    (1200.0_rk - 800.0_rk) * (1200.0_rk - co2_set)
            ENDIF

            ! Parameters for atmospheric CO2:
            ISMAXa    = 1.344_rk
            HEXPa     = 1.4614_rk
            CSTARa    = 585.0_rk

            ! For now, set CO2_Ci = 0.7d0 * CO2_Ca as recommended by Heald
            ! et al. (2009):
            CO2i      = 0.7_rk * co2_set

            ! Compute GAMMA_CO2:
            GAMMA_CO2 = ( ISMAXi -  ISMAXi * CO2i**HEXPi / &
                ( CSTARi**HEXPi + CO2i**HEXPi ) )  &
                * ( ISMAXa - ISMAXa * ( 0.7_rk * co2_set )**HEXPa / &
                ( CSTARa**HEXPa + ( 0.7_rk * co2_set )**HEXPa ) )

        ELSE

            ! No CO2 inhibition scheme is used; GAMMA_CO2 set to unity:
            GAMMA_CO2 = 1.0_rk

        ENDIF

    end function GET_GAMMA_CO2

    function GET_GAMMA_LEAFAGE(leafage_opt,LAIp,LAIc,tdays,Tt,Anew,Agro,Amat,Aold)  result( GAMMA_LEAFAGE )
        ! ROUTINE: GET_GAMMA_LEAFAGE
        !
        ! !DESCRIPTION: Function GET_GAMMA_LEAFAGE computes the leaf age activity factor
        !  associated with foliage fraction calculation.
        !
        !     leaf age response to Biogenic VOCs
        ! Revision: Sept 2023 Quazi Z. Rasool NOAA CSL/CIRES
        !----------------------------------------------------------------
        !
        !       GAMLA = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
        !       where Fnew = new foliage fraction
        !             Fgro = growing foliage fraction
        !             Fmat = mature foliage fraction
        !             Fold = old foliage fraction
        !             Anew = emission activity for new foliage
        !             Agro = emission activity for growing foliage
        !             Amat = emission activity for mature foliage
        !             Aold = emission activity for old foliage
        !           "Age class fractions are determined from LAI changes"
        !             LAIc = current Month's LAI (asuuming monthly LAI but can be
        !             customized as per tdays)
        !             LAIp = past Month's LAI
        !             tdays  = length of the time step (days) i.e. days in between LAIp and LAIc
        !             ti = days between budbreak and emission induction (calculated below)
        !             tm = days between budbreak and peak emission (Calculated below)
        !             Tt = 2-meter temperature (K) TEMP2 from the input model/obs
        !------------------------------------------------------------------------------
        ! !INTERFACE:

        ! !INPUT PARAMETERS:
        INTEGER,  INTENT(IN) :: leafage_opt       ! Option for leaf age emission factor calculation
        ! 0=On;
        ! 1 or  >1 =off i.e. GAMMA_LEAFAGE =1

        REAl(rk), INTENT(IN) :: tdays                      ! time step, number of days between Past and Current LAI inputs
        REAL(rk), INTENT(IN) :: Tt !AboveCanopy_TEMP(:)       ! Above canopy temperature (K), t2m or tmpsfc
        REAL(rk), INTENT(IN) :: LAIp         ! Past LAI [cm2/cm2]
        REAL(rk), INTENT(IN) :: LAIc         ! Current LAI [cm2/cm2]
        REAL(rk), INTENT(IN) :: Anew        ! Relative emiss factor (new leaves)
        REAL(rk), INTENT(IN) :: Agro        ! Relative emiss factor (growing leaves)
        REAL(rk), INTENT(IN) :: Amat        ! Relative emiss factor (mature leaves)
        REAL(rk), INTENT(IN) :: Aold        ! Relative emiss factor (old leaves)
        !
        ! !RETURN VALUE:
        REAL(rk) :: GAMMA_LEAFAGE             ! Leaf age activity factor [unitless]
        !
        ! !LOCAL VARIABLES:
        !INTEGER :: tdays                      ! time step
        REAL(rk) :: Fnew, Fgro                ! foliage fractions
        REAL(rk) :: Fmat, Fold
        REAL(rk) :: ti, tm
        !REAL(rk) :: Tt
        !
        ! !REMARKS:
        !  The function computes the leaf age activity factor based on BVOC's leaf age response.
        !
        ! !REVISION HISTORY:
        !  Adapted from HEMCO: https://github.com/geoschem/hemco for complete history
        !EOP
        !------------------------------------------------------------------------------
        !BOC

        !-----------------------
        ! Compute GAMMA_LEAFAGE
        !-----------------------


        !Tt = AboveCanopy_TEMP

        ! Calculate foliage fraction
        !-----------------------
        !Also, Compute ti and tm
        ! ti: number of days after budbreak required to induce emissions
        ! tm: number of days after budbreak required to reach peak emissions
        IF (LAIp < LAIc) THEN !(i.e. LAI has Increased)
            IF (Tt .le. 303.0_rk) THEN
                ti = 5.0_rk + 0.7_rk*(300.0_rk-Tt)
            ELSE
                ti = 2.9_rk
            ENDIF
            tm = 2.3_rk*ti
            !Fnew calculated
            IF (ti .ge. tdays) THEN
                Fnew = 1.0_rk - (LAIp/LAIc)
            ELSE
                Fnew = (ti/tdays) * ( 1.0_rk-(LAIp/LAIc) )
            ENDIF
            !Fmat calculated
            IF (tm .ge. tdays) THEN
                Fmat = LAIp/LAIc
            ELSE
                Fmat = (LAIp/LAIc) + ( (tdays-tm)/tdays ) * ( 1.0_rk-(LAIp/LAIc) )
            ENDIF

            Fgro = 1.0_rk - Fnew - Fmat
            Fold = 0.0_rk
        ELSEIF (LAIp == LAIc) THEN !If LAI remains same
            Fnew = 0.0_rk
            Fgro = 0.1_rk
            Fmat = 0.8_rk
            Fold = 0.1_rk
        ELSE
            Fnew = 0.0_rk
            Fgro = 0.0_rk
            Fold = ( LAIp-LAIc ) / LAIp
            Fmat = 1.0_rk-Fold
        ENDIF

        !-----------------------
        ! Compute GAMMA_AGE
        !-----------------------
        IF (leafage_opt .eq. 0) THEN
            ! Compute GAMMA_LEAFAGE
            GAMMA_LEAFAGE = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
            ! Prevent negative values
            GAMMA_LEAFAGE = MAX( GAMMA_LEAFAGE , 0.0_rk )
        ELSEIF (leafage_opt .eq. 1) THEN
            GAMMA_LEAFAGE = 1.0_rk
        ELSE
            GAMMA_LEAFAGE = 1.0_rk
        ENDIF

    end function GET_GAMMA_LEAFAGE




end module canopy_utils_mod
