program ode_logistic

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  real(wp), parameter :: y0 = 0.2_wp
  real(wp), parameter :: tol = 1.0e-8_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps, exact, err, c

  y = [y0]
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 2.0_wp
  h0 = 0.01_wp
  eps = 1.0e-10_wp

  call stiff3(n,fun,jac,x0,x1,h0,eps,w,y)

  c = (1.0_wp - y0)/y0
  exact = 1.0_wp/(1.0_wp + c*exp(-x1))
  err = abs(y(1) - exact)
  if (err > tol) then
    print '(A,ES12.4,A,ES12.4)', 'logistic error=', err, ' exceeds tolerance=', tol
    error stop 1
  end if

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = y(1)*(1.0_wp - y(1))
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = 1.0_wp - 2.0_wp*y(1)
  end subroutine

end program
