program ode_workspace

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: tol = 1.0e-12_wp
  real(wp) :: y_auto(n), y_work(n), w(n), x0, x1, h0, eps
  real(wp), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  integer :: stats_auto(3), stats_work(3)

  y_auto = [1.0_wp, -0.5_wp]
  y_work = y_auto
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.02_wp
  eps = 1.0e-8_wp

  allocate(rwork(n*(7 + 2*n)))
  allocate(iwork(n))

  call stiff3(n,fun,jac,x0,x1,h0,eps,w,y_auto,stats=stats_auto)

  h0 = 0.02_wp
  call stiff3(n,fun,jac,x0,x1,h0,eps,w,y_work,rwork,iwork,stats=stats_work)

  if (maxval(abs(y_auto - y_work)) > tol) then
    print '(A,2(1X,ES12.4))', 'solutions differ:', y_auto - y_work
    error stop 1
  end if

  if (any(stats_auto /= stats_work)) then
    print '(A,3(I0,1X),A,3(I0,1X))', &
      'stats mismatch auto [nfev njev nlu]: ', stats_auto, &
      ' work: ', stats_work
    error stop 1
  end if

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -2.0_wp*y(1) + y(2)
    f(2) = -3.0_wp*y(2)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -2.0_wp
    df(1,2) = 1.0_wp
    df(2,1) = 0.0_wp
    df(2,2) = -3.0_wp
  end subroutine

end program
