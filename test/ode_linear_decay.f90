program ode_linear_decay

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  real(wp), parameter :: lambda = 2.0_wp
  real(wp), parameter :: tol = 1.0e-8_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps, exact, err
  integer :: idid

  y = [1.0_wp]
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.01_wp
  eps = 1.0e-10_wp

  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,idid)
  if (idid /= 0) then
    print '(A,I0)', 'stiff3 failed with idid=', idid
    error stop 1
  end if

  exact = exp(-lambda*x1)
  err = abs(y(1) - exact)
  if (err > tol) then
    print '(A,ES12.4,A,ES12.4)', 'linear decay error=', err, ' exceeds tolerance=', tol
    error stop 1
  end if

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -lambda*y(1)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -lambda
  end subroutine

end program
