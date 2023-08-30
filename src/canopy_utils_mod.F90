module canopy_utils_mod

    use canopy_const_mod, ONLY: pi, rk, rearth    !constants for canopy models

    implicit none

    private
    public IntegrateTrapezoid,interp_linear1_internal,CalcPAI, &
        CalcDX,CalcFlameH,GET_GAMMA_CO2

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

    real(rk) function CalcFlameH(frp, dx, lu_opt, vtype) result( CalcFlameH_set )
        !! Approximates the Flame Height as a function of FRP intensity, grid cell distance (dx),
        !! and vegetation type (Based on Alexander and Cruz 2012).

        !!  Alexander, Martin E.; Cruz, Miguel G. 2012. Interdependencies between flame length and
        !!  fireline intensity in predicting crown fire initiation and crown scorch height.
        !!  International Journal of Wildland Fire 21(2):95-113.

        !! Assume Flame Length = Flame Height under calm winds

        real(rk), intent(in)  :: frp              !! Input Grid cell Fire Radiative Power (MW/cell)
        real(rk), intent(in)  :: dx               !! Input Grid cell length (m)
        integer,  intent(in)  :: lu_opt           !! Supported land use classifications
        integer,  intent(in)  :: vtype            !! Grid cell dominant vegetation type
        real(rk)              :: CalcFlameH_B59   !! Fireline intensity flame length (height) relationship (m): Byram (1959)
        real(rk)              :: CalcFlameH_A66_L !! Fireline intensity to flame height (m): Anderson et al. (1966)
        real(rk)              :: CalcFlameH_A66_D !! Fireline intensity to flame height (m): Anderson et al. (1966)
        real(rk)              :: CalcFlameH_N80   !! Fireline intensity to flame height (m): Nelson (1980)
        real(rk)              :: CalcFlameH_C83_H !! Fireline intensity to flame height (m): Clark (1983)
        real(rk)              :: CalcFlameH_C83_B !! Fireline intensity to flame height (m): Clark (1983)
        real(rk)              :: CalcFlameH_N86   !! Fireline intensity to flame height (m): Nelson and Adkins (1986)
        real(rk)              :: CalcFlameH_V86   !! Fireline intensity to flame height (m): van Wilgen (1986)
        real(rk)              :: CalcFlameH_B94   !! Fireline intensity to flame height (m): Burrows (1994)
        real(rk)              :: CalcFlameH_W96   !! Fireline intensity to flame height (m): Weise and Biging (1996)
        real(rk)              :: CalcFlameH_V98   !! Fireline intensity to flame height (m): Vega et al. (1998)
        real(rk)              :: CalcFlameH_C98   !! Fireline intensity to flame height (m): Catchpole et al. (1998)
        real(rk)              :: CalcFlameH_F00   !! Fireline intensity to flame height (m): Fernandes et al. (2000)
        real(rk)              :: CalcFlameH_B04   !! Fireline intensity to flame height (m): Butler et al. (2004)
        real(rk)              :: CalcFlameH_F09_H !! Fireline intensity to flame height (m): Fernandes et al. (2009)
        real(rk)              :: CalcFlameH_F09_B !! Fireline intensity to flame height (m): Fernandes et al. (2009)

        !!  -------------------------------------------------------------------------------------
        !!  Byram, GM (1959). Combustion of Forest Fuels. In Forest Fire: Control and Use.
        !!  (Ed. KP David) pp. 61-89.  McGraw Hill, New York, NY
        !Pine litter with grass understorey (Evergreen and Grasses)
        CalcFlameH_B59=0.0775_rk*((frp*1000.0_rk)/dx)**0.46_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Anderson HE, Brackebusch AP, Mutch RW, Rothermel RC (1966) Mechanisms of fire spread
        !!  research progress report 2. USDA Forest Service,
        !!  Intermountain Forest and Range Experiment Station, Research Paper
        !!  INT-28. (Ogden, UT)
        !Lodgepole pine slash (Evergreen)
        CalcFlameH_A66_L=0.074_rk*((frp*1000.0_rk)/dx)**0.651_rk
        !Douglas-fir slash  (Evergreen)
        CalcFlameH_A66_D=0.0447_rk*((frp*1000.0_rk)/dx)**0.67_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Nelson RM Jr (1980) Flame characteristics for fires in southern fuels. USDA
        !!  Forest Service, Southeastern Forest Experiment Station, Research Paper
        !!  SE-205. (Asheville, NC)
        !Southern USA Fuels (Mixed Forest)
        CalcFlameH_N80=0.0377_rk*((frp*1000.0_rk)/dx)**0.50_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Clark RG (1983) Threshold requirements for fire spread in grassland fuels.
        !!  PhD dissertation, Texas Tech University, Lubbock.
        !Grasslands (head fire) (Grasslands)
        CalcFlameH_C83_H=0.00015_rk*((frp*1000.0_rk)/dx)**1.75_rk
        !Grasslands (backfire) (Grasslands)
        CalcFlameH_C83_B=0.000722_rk*((frp*1000.0_rk)/dx)**0.99_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Nelson RM Jr, Adkins CW (1986) Flame characteristics of wind-driven
        !!  surface fires. Canadian Journal of Forest Research 16, 1293–1300.
        !!  doi:10.1139/X86-229
        !Litter and shrubs (Shrublands)
        CalcFlameH_N86=0.0475_rk*((frp*1000.0_rk)/dx)**0.493_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  van Wilgen BW (1986) A simple relationship for estimating the intensity of
        !!  fires in natural vegetation. South African Journal of Botany 52, 384–385
        !Fynbos shrublands (Shrublands)
        CalcFlameH_V86=0.0075_rk*((frp*1000.0_rk)/dx)**0.46_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Burrows ND (1994) Experimental development of a fire management model
        !!  for jarrah (Eucalyptus marginata Donn ex Sm.) forest. PhD thesis,
        !!  Australian National University, Canberra.
        !Eucalypt Forest (Evergreens)
        CalcFlameH_B94=0.0147_rk*((frp*1000.0_rk)/dx)**0.767_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Weise DR, Biging GS (1996) Effects of wind velocity and slope on flame
        !!  properties. Canadian Journal of Forest Research 26, 1849–1858.
        !!  doi:10.1139/X26-210
        !Excelsior (Deciduous)
        CalcFlameH_W96=0.016_rk*((frp*1000.0_rk)/dx)**0.7_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Vega JA, Cuinas P, Fonturbel T, Perez-Gorostiaga P, Fernandez C (1998)
        !!  Predicting fire behaviour in Galician (NW Spain) shrubland fuel
        !!  complexes. In ‘Proceedings of 3rd International Conference on Forest
        !!  Flame length and fireline intensity interdependences.
        !Shrublands
        CalcFlameH_V98=0.087_rk*((frp*1000.0_rk)/dx)**0.493_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  CatchpoleWR, Bradstock RA, Choate J, Fogarty LG, Gellie N, McCarthy G,
        !!  McCaw WL, Marsden-Smedley JB, Pearce G (1998) Cooperative
        !!  development of equations for heathland fire behaviour. In ‘Proceedings
        !!  of 3rd International Conference on Forest Fire Research and 14th
        !!  Conference on Fire and Forest Meteorology, Volume II’, 16–20
        !!  November 1998, Luso–Coimbra, Portugal. (Ed. DX Viegas)
        !!  pp. 631–645. (University of Coimbra: Coimbra, Portugal)
        !Shrublands
        CalcFlameH_C98=0.0325_rk*((frp*1000.0_rk)/dx)**0.56_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Fernandes PM, Catchpole WR, Rego FC (2000) Shrubland fire behaviour
        !!  modelling with microplot data.Canadian Journal of Forest Research 30,
        !!  889–899. doi:10.1139/X00-012
        !Shrublands
        CalcFlameH_F00=0.0516_rk*((frp*1000.0_rk)/dx)**0.453_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Butler BW, Finney MA, Andrews PL, Albini FA (2004) A radiation-driven
        !!  model of crown fire spread. Canadian Journal of Forest Research 34,
        !!  1588–1599. doi:10.1139/X04-074
        !Jack Pine Forest - Crown Fire (Evergreen)
        CalcFlameH_B04=0.0175_rk*((frp*1000.0_rk)/dx)**0.667_rk
        !!  -------------------------------------------------------------------------------------

        !!  -------------------------------------------------------------------------------------
        !!  Fernandes PM, Botelho HS, Rego FC, Loureiro C (2009) Empirical
        !!  modelling of surface fire behaviour in maritime pine stands.
        !!  International Journal of Wildland Fire 18, 698–710. doi:10.1071/
        !!  WF08023
        !Maritime Pine (head fire) (Evergreen)
        CalcFlameH_F09_H=0.045_rk*((frp*1000.0_rk)/dx)**0.543_rk
        !Maritime Pine (backfire) (Evergreen)
        CalcFlameH_F09_B=0.029_rk*((frp*1000.0_rk)/dx)**0.724_rk
        !!  -------------------------------------------------------------------------------------

        if (lu_opt .eq. 0 .or. lu_opt .eq. 1) then !VIIRS or MODIS LU types
            if (vtype .ge. 1 .and. vtype .le. 2) then !VIIRS/MODIS Cat 1-2/Evergreen Needleleaf & Broadleaf
                CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_A66_L+CalcFlameH_A66_D+CalcFlameH_B94+ &
                    CalcFlameH_F09_H+CalcFlameH_F09_B) / 6.0_rk
                if (((frp*1000.0_rk)/dx) .ge. 1700.0_rk ) then !Evergreen forest crowning likely
                    CalcFlameH_set=CalcFlameH_B04
                end if
            else if (vtype .ge. 3 .and. vtype .le. 4) then !VIIRS/MODIS Cat 3-4 Deciduous Needleleaf and  Broadleaf
                CalcFlameH_set=(CalcFlameH_N80+CalcFlameH_W96) / 2.0_rk
            else if (vtype .eq. 5) then !VIIRS/MODIS Cat 5 Mixed Forests
                CalcFlameH_set=CalcFlameH_N80
            else if (vtype .ge. 6 .and. vtype .le. 7) then !VIIRS/MODIS Cat 6-7 Shrublands
                CalcFlameH_set=(CalcFlameH_N86+CalcFlameH_V86+CalcFlameH_V98+CalcFlameH_C98+ &
                    CalcFlameH_F00) / 5.0_rk
            else if (vtype .ge. 8 .and. vtype .le. 10) then !VIIRS/MODIS Cat 8-10 Savannas and Grasslands
                CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_C83_H+CalcFlameH_C83_B) / 3.0_rk
            else if (vtype .ge. 12) then !VIIRS/MODIS Cat 12 Croplands
                CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_C83_H+CalcFlameH_C83_B) / 3.0_rk
            else
                CalcFlameH_set=(CalcFlameH_B59+CalcFlameH_A66_L+CalcFlameH_A66_D+CalcFlameH_N80+ &
                    CalcFlameH_C83_H+CalcFlameH_C83_B+CalcFlameH_N86+CalcFlameH_V86+ &
                    CalcFlameH_B94+CalcFlameH_W96+CalcFlameH_V98+CalcFlameH_C98+ &
                    CalcFlameH_F00+CalcFlameH_B04+CalcFlameH_F09_H+ &
                    CalcFlameH_F09_B) / 16.0_rk
            end if
        else
            write(*,*)  'Wrong LU_OPT choice of ', LU_OPT, 'in namelist, only VIIRS/MODIS available right now...exiting'
            call exit(2)
        end if
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

end module canopy_utils_mod
