program ode_hmax

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  real(wp), parameter :: tol = 1.0e-10_wp
  real(wp), parameter :: tol_x = 1.0e-12_wp
  real(wp), parameter :: hmax_cap = 0.1_wp
  real(wp), parameter :: h0_start = 0.5_wp
  real(wp) :: y_default(n), y_zero(n), y_capped(n), y_stop(n), w(n), x0, x1, h0, eps_test, x
  real(wp) :: max_h_seen, x_last_accepted
  integer :: stats_default(6), stats_zero(6), stats_capped(6), stats_stop(6)

  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  eps_test = 1.0e-4_wp

  y_default = [1.0_wp]
  h0 = h0_start
  x = x0
  call stiff3(n,fun,x,y_default,x1,jac,h0,eps_test,w,stats=stats_default)
  if (abs(x - x1) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'default run x mismatch: ', x, ' expected ', x1
    error stop 1
  end if

  y_zero = [1.0_wp]
  h0 = h0_start
  x = x0
  call stiff3(n,fun,x,y_zero,x1,jac,h0,eps_test,w,stats=stats_zero,hmax=0.0_wp)
  if (abs(x - x1) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'hmax=0 run x mismatch: ', x, ' expected ', x1
    error stop 1
  end if

  if (any(stats_zero /= stats_default)) then
    print '(A,6(I0,1X),A,6(I0,1X))', &
      'expected hmax=0 stats [nacc nrej nfev njev nlu nsol] to match default. got ', stats_zero, &
      ' expected ', stats_default
    error stop 1
  end if

  if (abs(y_zero(1) - y_default(1)) > tol) then
    print '(A,ES12.4)', 'expected hmax=0 solution to match default, error=', abs(y_zero(1) - y_default(1))
    error stop 1
  end if

  y_capped = [1.0_wp]
  h0 = h0_start
  max_h_seen = 0.0_wp
  x = x0
  call stiff3(n,fun,x,y_capped,x1,jac,h0,eps_test,w,solout=out_cap,stats=stats_capped,hmax=hmax_cap)
  if (abs(x - x1) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'capped run x mismatch: ', x, ' expected ', x1
    error stop 1
  end if

  if (max_h_seen > hmax_cap + tol) then
    print '(A,ES12.4,A,ES12.4)', 'observed half-step exceeds hmax cap: ', max_h_seen, ' > ', hmax_cap
    error stop 1
  end if

  y_stop = [1.0_wp]
  h0 = h0_start
  x_last_accepted = x0
  x = x0
  call stiff3(n,fun,x,y_stop,x1,jac,h0,eps_test,w,solout=out_stop,stats=stats_stop,hmax=hmax_cap)

  if (abs(x - x_last_accepted) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'interrupt run x not last accepted step: ', x, ' expected ', x_last_accepted
    error stop 1
  end if

  if (x <= x0 + tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'interrupt run x too close to start: ', x, ' start ', x0
    error stop 1
  end if

  if (x >= x1 - tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'interrupt run x too close to end: ', x, ' end ', x1
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

  subroutine out_cap(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn
    max_h_seen = max(max_h_seen,abs(x - xold)/2.0_wp)
  end subroutine

  subroutine out_stop(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn
    x_last_accepted = x
    if (nr >= 3) irtrn = -1
  end subroutine

end program
