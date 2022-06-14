module Canopy_Utils

implicit none

private
public IntegrateTrapezoid

contains 

function IntegrateTrapezoid(x, y)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real, intent(in)  :: x(:)                !! Variable x
    real, intent(in)  :: y(size(x))          !! Function y(x)
    real              :: IntegrateTrapezoid  !! Integral ∫y(x)·dx
! Integrate using the trapezoidal rule
    associate(n => size(x))
      IntegrateTrapezoid = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
end function

end module Canopy_Utils
