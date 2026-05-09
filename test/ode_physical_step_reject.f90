program ode_physical_step_reject

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  real(wp), parameter :: y0 = 1.0_wp
  real(wp), parameter :: tol = 1.0e-7_wp
  real(wp) :: y_ref(n), y_retry(n), y_fail(n), w(n), x0, x1, h0, eps
  integer :: stats_ref(3), stats_retry(3), stats_fail(3)
  integer :: reject_once_count, reject_fail_count
  logical :: rejected_once

  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 0.5_wp
  eps = 1.0e-7_wp

  y_ref = [y0]
  h0 = 0.1_wp
  call stiff3(n,fun,x0,y_ref,x1,jac,h0,eps,w,stats=stats_ref)

  y_retry = [y0]
  h0 = 0.1_wp
  reject_once_count = 0
  rejected_once = .false.
  call stiff3(n,fun,x0,y_retry,x1,jac,h0,eps,w,solout=out_reject_once,stats=stats_retry)

  if (.not. rejected_once) then
    print '(A)', 'expected one callback-driven step rejection'
    error stop 1
  end if
  if (reject_once_count < 2) then
    print '(A,I0)', 'expected callback to be re-entered after rejection, count=', reject_once_count
    error stop 1
  end if
  if (abs(y_retry(1) - y_ref(1)) > tol) then
    print '(A,ES12.4)', 'retry solution drift exceeds tolerance=', abs(y_retry(1) - y_ref(1))
    error stop 1
  end if
  if (stats_retry(1) <= stats_ref(1)) then
    print '(A,2(I0,1X))', 'expected extra rhs work after retry, got nfev retry/ref=', stats_retry(1), stats_ref(1)
    error stop 1
  end if

  y_fail = [y0]
  h0 = 0.1_wp
  reject_fail_count = 0
  call stiff3(n,fun,x0,y_fail,x1,jac,h0,eps,w,solout=out_reject_always,stats=stats_fail)

  if (reject_fail_count < 2) then
    print '(A,I0)', 'expected repeated rejections before stop, count=', reject_fail_count
    error stop 1
  end if
  if (abs(y_fail(1) - y0) > tol) then
    print '(A,ES12.4)', 'expected rollback to initial state after failed retries, error=', abs(y_fail(1) - y0)
    error stop 1
  end if
  if (stats_fail(1) == 0) then
    print '(A)', 'expected solver work stats before retry stop'
    error stop 1
  end if

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -10.0_wp*y(1)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -10.0_wp
  end subroutine

  subroutine out_reject_once(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn
    reject_once_count = reject_once_count + 1
    if (.not. rejected_once) then
      irtrn = -1
      rejected_once = .true.
    end if
  end subroutine

  subroutine out_reject_always(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn
    reject_fail_count = reject_fail_count + 1
    irtrn = -1
  end subroutine

end program
