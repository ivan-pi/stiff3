program ode_dense_output_integral

  use stiff3_solver, only: stiff3, stiff3_interp, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: x0 = 0.0_wp
  real(wp), parameter :: x1 = 1.0_wp
  real(wp), parameter :: hmax = 0.05_wp
  real(wp), parameter :: dx_fixed = 1.0e-3_wp
  real(wp), parameter :: tol_match = 2.0e-4_wp
  real(wp), parameter :: tol_exact = 2.0e-4_wp
  real(wp), parameter :: tol_x = 1.0e-12_wp
  real(wp) :: y(n), w(n), h0, eps, x
  real(wp) :: rwork(n*(7 + 2*n))
  integer :: iwork(n), idid
  integer :: nsteps_checked
  real(wp) :: adaptive_integral, interp_integral, exact_integral
  real(wp) :: xprev, yprev1, xfixed_prev, yfixed_prev, xfixed_next

  y = [1.0_wp, 1.0_wp]
  w = 1.0_wp
  h0 = 0.05_wp
  eps = 1.0e-8_wp
  nsteps_checked = 0
  adaptive_integral = 0.0_wp
  interp_integral = 0.0_wp
  exact_integral = 1.0_wp - exp(-x1)
  xprev = x0
  yprev1 = y(1)
  xfixed_prev = x0
  yfixed_prev = y(1)
  xfixed_next = x0 + dx_fixed

  x = x0
  call stiff3(n,fun,x,y,x1,jac,h0,eps,w,rwork,iwork,idid,solout=accumulate_integrals,hmax=hmax)
  if (idid /= 0) then
    print '(A,I0)', 'expected successful solve, got idid=', idid
    error stop 1
  end if

  if (nsteps_checked < 2) then
    print '(A,I0)', 'expected integral test to exercise multiple steps, got ', nsteps_checked
    error stop 1
  end if

  if (abs(xfixed_prev - x1) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'fixed-grid integration stopped at ', xfixed_prev, ' expected ', x1
    error stop 1
  end if

  if (abs(x - x1) > tol_x) then
    print '(A,ES12.4,A,ES12.4)', 'solver x did not end at x1: ', x, ' expected ', x1
    error stop 1
  end if

  if (abs(interp_integral - adaptive_integral) > tol_match) then
    print '(A,ES12.4,A,ES12.4)', 'dense-output integral mismatch: ', interp_integral, ' vs ', adaptive_integral
    error stop 1
  end if

  if (abs(interp_integral - exact_integral) > tol_exact) then
    print '(A,ES12.4,A,ES12.4)', 'dense-output integral inaccurate: ', interp_integral, ' exact ', exact_integral
    error stop 1
  end if

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -y(1)
    f(2) = -50.0_wp*y(2)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -1.0_wp
    df(1,2) = 0.0_wp
    df(2,1) = 0.0_wp
    df(2,2) = -50.0_wp
  end subroutine

  subroutine accumulate_integrals(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    ! This callback compares two integral estimates of y(1) over [x0, x1].
    ! It carries parent-scope state through xprev/yprev1 for the native
    ! accepted-step trapezoid rule and through xfixed_prev/yfixed_prev/
    ! xfixed_next for the fine fixed-grid trapezoid rule sampled with
    ! stiff3_interp. The explicit workspace rwork is also inherited
    ! deliberately so stiff3_interp can access the accepted-step data.
    real(wp) :: yinterp

    if (abs(xold - xprev) > tol_x) then
      print '(A,I0,A,ES12.4,A,ES12.4)', &
        'accepted-step x mismatch on step ', nr, '. got=', xold, ' expected=', xprev
      error stop 1
    end if

    adaptive_integral = adaptive_integral + 0.5_wp*(x - xold)*(yprev1 + y(1))

    do while (xfixed_next <= x + tol_x)
      call stiff3_interp(xold,x,y,rwork,xfixed_next,1,yinterp)
      interp_integral = interp_integral + 0.5_wp*(xfixed_next - xfixed_prev)*(yfixed_prev + yinterp)
      xfixed_prev = xfixed_next
      yfixed_prev = yinterp
      xfixed_next = xfixed_next + dx_fixed
    end do

    nsteps_checked = nsteps_checked + 1
    xprev = x
    yprev1 = y(1)
  end subroutine

end program
