program ode_idid

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  real(wp), parameter :: rosenbrock_a = 0.4358665215084589_wp
  real(wp), parameter :: x0_underflow = 1.0e12_wp
  real(wp), parameter :: dx_underflow = 8.0_wp*spacing(x0_underflow)
  real(wp), parameter :: omega = 1.0e10_wp

  call check_lu_failure()
  call check_bisection_underflow()

contains

  subroutine check_lu_failure()
    real(wp) :: y(1), w(1), x, x1, h0, eps
    integer :: idid

    y = [1.0_wp]
    w = [1.0_wp]
    x = 0.0_wp
    x1 = 0.2_wp
    h0 = 0.1_wp
    eps = 1.0e-6_wp

    call stiff3(1,fun_lu,x,y,x1,jac_lu,h0,eps,w,idid)

    if (idid /= -1) then
      print '(A,I0)', 'expected LU failure idid=-1, got ', idid
      error stop 1
    end if

    if (x /= 0.0_wp) then
      print '(A,ES12.4)', 'expected LU failure to leave x unchanged, got ', x
      error stop 1
    end if
  end subroutine

  subroutine check_bisection_underflow()
    real(wp) :: y(2), w(2), x, x1, h0, eps
    integer :: idid

    y = [1.0_wp, 0.0_wp]
    w = [1.0_wp, 1.0_wp]
    x = x0_underflow
    x1 = x0_underflow + dx_underflow
    h0 = 0.5_wp*dx_underflow
    eps = 1.0e-14_wp

    call stiff3(2,fun_underflow,x,y,x1,jac_underflow,h0,eps,w,idid)

    if (idid /= -3) then
      print '(A,I0)', 'expected step-size underflow idid=-3, got ', idid
      error stop 1
    end if

    if (abs(x - x0_underflow) > spacing(x0_underflow)) then
      print '(A,ES12.4,A,ES12.4)', 'expected underflow to leave x near start: ', x, ' start ', x0_underflow
      error stop 1
    end if
  end subroutine

  subroutine fun_lu(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f = 0.0_wp
  end subroutine

  subroutine jac_lu(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = 1.0_wp/(0.2_wp*rosenbrock_a)
  end subroutine

  subroutine fun_underflow(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = omega*y(2)
    f(2) = -omega*y(1)
  end subroutine

  subroutine jac_underflow(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = 0.0_wp
    df(1,2) = omega
    df(2,1) = -omega
    df(2,2) = 0.0_wp
  end subroutine

end program
