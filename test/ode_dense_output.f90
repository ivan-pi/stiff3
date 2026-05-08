program ode_dense_output

  use stiff3_solver, only: stiff3, stiff3_interp_component, stiff3_interp_vector, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: tol_mid = 1.0e-6_wp
  real(wp), parameter :: tol_end = 1.0e-12_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  real(wp), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  logical :: checked_dense_output

  y = [1.0_wp, 1.0_wp]
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.5_wp
  eps = 1.0e-8_wp
  checked_dense_output = .false.

  allocate(rwork(n*(8 + 2*n)))
  allocate(iwork(n))

  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,rwork,iwork,solout=check_dense_output)

  if (.not. checked_dense_output) then
    print '(A)', 'dense-output callback was not exercised'
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

  subroutine check_dense_output(nr,xold,x,y,iha,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    real(wp) :: xmid, yscalar
    real(wp) :: ymid(n), ystart(n), yend(n)
    real(wp) :: exact_mid(n)

    if (nr /= 1) return

    if (iha < 1) then
      print '(A,I0)', 'expected first accepted step to include at least one bisection, got iha=', iha
      error stop 1
    end if

    xmid = 0.5_wp*(xold + x)
    exact_mid = exact_state(xmid)

    call stiff3_interp_vector(xold,x,y,rwork,xmid,ymid)
    call stiff3_interp_component(xold,x,y,rwork,xmid,2,yscalar)

    if (maxval(abs(ymid - exact_mid)) > tol_mid) then
      print '(A,2(1X,ES12.4),A,2(1X,ES12.4))', &
        'midpoint interpolation mismatch. got=', ymid, ' expected=', exact_mid
      error stop 1
    end if

    if (abs(yscalar - ymid(2)) > tol_end) then
      print '(A,ES12.4,A,ES12.4)', 'scalar/vector interpolation mismatch: ', yscalar, ' vs ', ymid(2)
      error stop 1
    end if

    call stiff3_interp_vector(xold,x,y,rwork,xold,ystart)
    if (maxval(abs(ystart - [1.0_wp, 1.0_wp])) > tol_end) then
      print '(A,2(1X,ES12.4),A)', 'step-start interpolation mismatch: ', ystart, ' expected initial state'
      error stop 1
    end if

    call stiff3_interp_vector(xold,x,y,rwork,x,yend)
    if (maxval(abs(yend - y)) > tol_end) then
      print '(A,2(1X,ES12.4),A,2(1X,ES12.4))', 'step-end interpolation mismatch. got=', yend, ' expected=', y
      error stop 1
    end if

    checked_dense_output = .true.
    irtrn = -1
  end subroutine

  function exact_state(x) result(yexact)
    real(wp), intent(in) :: x
    real(wp) :: yexact(n)

    yexact(1) = exp(-x)
    yexact(2) = exp(-50.0_wp*x)
  end function

end program
