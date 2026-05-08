program ode_stats

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  real(wp), parameter :: lambda = 2.0_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  integer :: stats(3)

  y = [1.0_wp]
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.01_wp
  eps = 1.0e-10_wp

  call stiff3(n,fun,jac,x0,x1,h0,eps,w,y,stats=stats)

  if (any(stats <= 0)) then
    print '(A,3(I0,1X))', 'expected positive stats [nfev njev nlu], got: ', stats
    error stop 1
  end if

  if (stats(1) /= stats(2) + stats(3)) then
    print '(A,3(I0,1X))', 'inconsistent stats [nfev njev nlu], got: ', stats
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
