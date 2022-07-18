module canopy_utils_mod

use canopy_const_mod, ONLY: pi, rk    !constants for canopy models
        
implicit none

private
public IntegrateTrapezoid,interp_linear1_internal,CalcPAI

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


        real(rk), intent(in)  :: fch                 !! Input Grid cell canopy height (m
        real(rk), intent(in)  :: ffrac               !! Input Grid cell forest fraction
        real(rk)              :: CalcPAI             !! Calculated Plant area index (PAI)

        CalcPAI=( (fch*(ffrac/3.0_rk)*10.6955_rk) / (2.0_rk * pi) ) * ffrac !Massman PAI calculation (Eq. 19)
        
      end function 

end module canopy_utils_mod
