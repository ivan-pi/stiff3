program ode_workspace

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp) :: y_auto(n), y_work(n), w(n), x0, x1, h0, eps, x
  real(wp), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  integer :: stats_auto(6), stats_work(6), idid_auto, idid_work

  y_auto = [1.0_wp, -0.5_wp]
  y_work = y_auto
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.02_wp
  eps = 1.0e-8_wp

  allocate(rwork(n*(7 + 2*n)))
  allocate(iwork(n))

  x = x0
  call stiff3(n,fun,x,y_auto,x1,jac,h0,eps,w,idid_auto,stats=stats_auto)
  if (idid_auto /= 0) then
    print '(A,I0)', 'auto workspace expected idid=0, got ', idid_auto
    error stop 1
  end if
  if (x /= x1) then
    print '(A,ES12.4,A,ES12.4)', 'auto workspace x mismatch: ', x, ' expected ', x1
    error stop 1
  end if

  h0 = 0.02_wp
  x = x0
  call stiff3(n,fun,x,y_work,x1,jac,h0,eps,w,rwork,iwork,idid_work,stats=stats_work)
  if (idid_work /= 0) then
    print '(A,I0)', 'explicit workspace expected idid=0, got ', idid_work
    error stop 1
  end if
  if (x /= x1) then
    print '(A,ES12.4,A,ES12.4)', 'explicit workspace x mismatch: ', x, ' expected ', x1
    error stop 1
  end if

  if (any(y_auto /= y_work)) then
    print '(A,2(1X,ES12.4),A,2(1X,ES12.4))', &
      'solutions differ. auto=', y_auto, ' work=', y_work
    error stop 1
  end if

  if (any(stats_auto /= stats_work)) then
    print '(A,6(I0,1X),A,6(I0,1X))', &
      'stats mismatch auto [nacc nrej nfev njev nlu nsol]: ', stats_auto, &
      ' work: ', stats_work
    error stop 1
  end if

  if (idid_auto /= idid_work) then
    print '(A,2(I0,1X))', 'idid mismatch auto/work: ', idid_auto, idid_work
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
