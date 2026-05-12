program ode_hmax

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  real(wp), parameter :: tol = 1.0e-10_wp
  real(wp), parameter :: tol_x = 1.0e-12_wp
  real(wp), parameter :: hmax_cap = 0.1_wp
  real(wp), parameter :: h0_start = 0.5_wp
  real(wp) :: y_default(n), y_zero(n), y_capped(n), y_stop(n), y_rhs_stop(n), w(n), x0, x1, h0, eps_test, x
  real(wp) :: max_h_seen, x_last_accepted, x_last_rhs_accepted, y_last_rhs_accepted(n), rhs_stop_threshold
  integer :: stats_default(6), stats_zero(6), stats_capped(6), stats_stop(6), stats_rhs_stop(6)
  integer :: idid_default, idid_zero, idid_capped, idid_stop, idid_rhs_stop
  logical :: enable_rhs_interrupt

  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  eps_test = 1.0e-4_wp
  enable_rhs_interrupt = .false.

  y_default = [1.0_wp]
  h0 = h0_start
  x = x0
  call stiff3(n,fun,x,y_default,x1,jac,h0,eps_test,w,idid_default,stats=stats_default)
  if (idid_default /= 0) then
    print '(A,I0)', 'default run expected idid=0, got ', idid_default
    error stop 1
  end if
  if (abs(x - x1) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'default run x mismatch: ', x, ' expected ', x1
    error stop 1
  end if

  y_zero = [1.0_wp]
  h0 = h0_start
  x = x0
  call stiff3(n,fun,x,y_zero,x1,jac,h0,eps_test,w,idid_zero,stats=stats_zero,hmax=0.0_wp)
  if (idid_zero /= 0) then
    print '(A,I0)', 'hmax=0 run expected idid=0, got ', idid_zero
    error stop 1
  end if
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
  call stiff3(n,fun,x,y_capped,x1,jac,h0,eps_test,w,idid_capped,solout=out_cap,stats=stats_capped,hmax=hmax_cap)
  if (idid_capped /= 0) then
    print '(A,I0)', 'capped run expected idid=0, got ', idid_capped
    error stop 1
  end if
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
  call stiff3(n,fun,x,y_stop,x1,jac,h0,eps_test,w,idid_stop,solout=out_stop,stats=stats_stop,hmax=hmax_cap)
  if (idid_stop /= -2) then
    print '(A,I0)', 'interrupt run expected idid=-2, got ', idid_stop
    error stop 1
  end if

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

  y_rhs_stop = [1.0_wp]
  h0 = h0_start
  x_last_rhs_accepted = x0
  y_last_rhs_accepted = y_rhs_stop
  ! Chosen so the first accepted step completes, but the next step triggers
  ! inside an RHS evaluation and should therefore return the prior accepted state.
  rhs_stop_threshold = 0.75_wp
  enable_rhs_interrupt = .true.
  x = x0
  call stiff3(n,fun,x,y_rhs_stop,x1,jac,h0,eps_test,w,idid_rhs_stop,solout=out_rhs_track,stats=stats_rhs_stop,hmax=hmax_cap)
  enable_rhs_interrupt = .false.
  if (idid_rhs_stop /= -2) then
    print '(A,I0)', 'rhs interrupt run expected idid=-2, got ', idid_rhs_stop
    error stop 1
  end if
  if (abs(x - x_last_rhs_accepted) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'rhs interrupt run x not last accepted step: ', x, ' expected ', x_last_rhs_accepted
    error stop 1
  end if
  if (abs(y_rhs_stop(1) - y_last_rhs_accepted(1)) > tol) then
    print '(A,ES12.4,A,ES12.4)', &
      'rhs interrupt run y not last accepted value: ', y_rhs_stop(1), ' expected ', y_last_rhs_accepted(1)
    error stop 1
  end if
  if (x <= x0 + tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'rhs interrupt run x too close to start: ', x, ' start ', x0
    error stop 1
  end if

contains

  subroutine fun(n,y,f, ires)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    integer, intent(inout) :: ires
    f(1) = -y(1)
    if (enable_rhs_interrupt .and. y(1) < rhs_stop_threshold) ires = -1
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

  subroutine out_rhs_track(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn
    x_last_rhs_accepted = x
    y_last_rhs_accepted = y
  end subroutine

end program
