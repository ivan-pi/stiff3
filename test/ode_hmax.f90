program ode_hmax

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  real(wp), parameter :: tol = 1.0e-10_wp
  real(wp) :: y_default(n), y_zero(n), y_capped(n), w(n), x0, x1, h0, eps
  integer :: stats_default(3), stats_zero(3), stats_capped(3)

  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.5_wp
  eps = 1.0e6_wp

  y_default = [1.0_wp]
  call stiff3(n,fun,jac,x0,x1,h0,eps,w,y_default,stats=stats_default)

  y_zero = [1.0_wp]
  call stiff3(n,fun,jac,x0,x1,h0,eps,w,y_zero,stats=stats_zero,hmax=0.0_wp)

  if (any(stats_zero /= stats_default)) then
    print '(A,3(I0,1X),A,3(I0,1X))', &
      'expected hmax=0 stats [nfev njev nlu] to match default. got ', stats_zero, &
      ' expected ', stats_default
    error stop 1
  end if

  if (abs(y_zero(1) - y_default(1)) > tol) then
    print '(A,ES12.4)', 'expected hmax=0 solution to match default, error=', abs(y_zero(1) - y_default(1))
    error stop 1
  end if

  y_capped = [1.0_wp]
  call stiff3(n,fun,jac,x0,x1,h0,eps,w,y_capped,stats=stats_capped,hmax=0.1_wp)

  if (stats_capped(3) <= stats_default(3)) then
    print '(A,3(I0,1X),A,3(I0,1X))', &
      'expected hmax cap to increase LU count [nfev njev nlu]. capped=', stats_capped, &
      ' default=', stats_default
    error stop 1
  end if

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -y(1)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -1.0_wp
  end subroutine

end program
