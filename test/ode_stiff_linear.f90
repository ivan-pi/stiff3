program ode_stiff_linear

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: tol = 2.0e-7_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  real(wp) :: exact1, exact2, err

  y = [2.0_wp, 1.0_wp]
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 0.1_wp
  h0 = 1.0e-3_wp
  eps = 1.0e-10_wp

  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w)

  exact2 = exp(x1)
  exact1 = (1000.0_wp/1001.0_wp)*exact2 + (1002.0_wp/1001.0_wp)*exp(-1000.0_wp*x1)
  err = max(abs(y(1) - exact1),abs(y(2) - exact2))
  if (err > tol) then
    print '(A,ES12.4,A,ES12.4)', 'stiff linear error=', err, ' exceeds tolerance=', tol
    error stop 1
  end if

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -1000.0_wp*(y(1) - y(2))
    f(2) = y(2)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -1000.0_wp
    df(1,2) = 1000.0_wp
    df(2,1) = 0.0_wp
    df(2,2) = 1.0_wp
  end subroutine

end program
