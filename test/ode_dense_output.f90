program ode_dense_output

  use stiff3_solver, only: stiff3, stiff3_interp, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: tol_mid = 1.0e-6_wp
  real(wp), parameter :: tol_end = 1.0e-12_wp
  real(wp), parameter :: hmax = 0.05_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  real(wp) :: rwork(n*(7 + 2*n))
  integer :: iwork(n), idid
  integer :: nsteps_checked
  real(wp) :: yprev(n), xprev

  y = [1.0_wp, 1.0_wp]
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.5_wp
  eps = 1.0e-8_wp
  nsteps_checked = 0
  xprev = x0
  yprev = y

  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,rwork,iwork,idid,solout=check_dense_output,hmax=hmax)
  if (idid /= 0) then
    print '(A,I0)', 'expected successful solve, got idid=', idid
    error stop 1
  end if

  if (nsteps_checked < 2) then
    print '(A,I0)', 'expected dense-output test to exercise multiple steps, got ', nsteps_checked
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

  function exact_state(x) result(yexact)
    real(wp), intent(in) :: x
    real(wp) :: yexact(n)

    yexact(1) = exp(-x)
    yexact(2) = exp(-50.0_wp*x)
  end function

  subroutine check_dense_output(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    ! This callback checks dense output on every accepted step over [xold, x].
    ! It carries parent-scope state through xprev and yprev, which hold the
    ! previous accepted endpoint so the next callback can verify continuity at
    ! the interpolant start. The explicit workspace rwork is also inherited
    ! from the parent scope on purpose so stiff3_interp can read the accepted
    ! step data stored by stiff3_work.
    real(wp) :: xmid, yscalar
    real(wp) :: ymid(n), ystart(n), yend(n)
    real(wp) :: exact_mid(n), expected_start(n)

    xmid = 0.5_wp*(xold + x)
    exact_mid = exact_state(xmid)

    call stiff3_interp(xold,x,y,rwork,xmid,ymid)
    call stiff3_interp(xold,x,y,rwork,xmid,2,yscalar)

    if (maxval(abs(ymid - exact_mid)) > tol_mid) then
      print '(A,2(1X,ES12.4),A,2(1X,ES12.4))', &
        'midpoint interpolation mismatch. got=', ymid, ' expected=', exact_mid
      error stop 1
    end if

    if (abs(yscalar - ymid(2)) > tol_end) then
      print '(A,ES12.4,A,ES12.4)', 'scalar/vector interpolation mismatch: ', yscalar, ' vs ', ymid(2)
      error stop 1
    end if

    if (nr == 1) then
      expected_start = exact_state(xold)
    else
      expected_start = yprev
      if (abs(xold - xprev) > tol_end) then
        print '(A,I0,A,ES12.4,A,ES12.4)', &
          'step-start x-coordinate mismatch on step ', nr, '. got=', xold, ' expected=', xprev
        error stop 1
      end if
    end if

    call stiff3_interp(xold,x,y,rwork,xold,ystart)
    if (maxval(abs(ystart - expected_start)) > tol_end) then
      print '(A,I0,A,2(1X,ES12.4),A,2(1X,ES12.4))', &
        'step-start interpolation mismatch on step ', nr, '. got=', ystart, ' expected=', expected_start
      error stop 1
    end if

    call stiff3_interp(xold,x,y,rwork,x,yend)
    if (maxval(abs(yend - y)) > tol_end) then
      print '(A,2(1X,ES12.4),A,2(1X,ES12.4))', 'step-end interpolation mismatch. got=', yend, ' expected=', y
      error stop 1
    end if

    nsteps_checked = nsteps_checked + 1
    xprev = x
    yprev = y
  end subroutine

end program
