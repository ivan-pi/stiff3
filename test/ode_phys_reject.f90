program ode_phys_reject

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1
  integer, parameter :: reject_injected = 3
  real(wp), parameter :: tol_x = 1.0e-12_wp
  real(wp) :: y(n), w(n), x, x0, x1, h0, eps
  integer :: idid, stats(6)
  logical :: always_reject
  integer :: reject_remaining

  call check_recovery_path()
  call check_reject_limit_path()

contains

  subroutine check_recovery_path()
    x0 = 0.0_wp
    x1 = 1.0_wp
    x = x0
    y = [1.0_wp]
    w = [1.0_wp]
    h0 = 0.25_wp
    eps = 1.0e-6_wp

    always_reject = .false.
    reject_remaining = reject_injected
    call stiff3(n,fun,x,y,x1,jac,h0,eps,w,idid,stats=stats)

    if (idid /= 0) then
      print '(A,I0)', 'recovery case expected idid=0, got ', idid
      error stop 1
    end if
    if (abs(x - x1) > tol_x) then
      print '(A,ES12.4,A,ES12.4)', 'recovery case expected x=x1, got ', x, ' expected ', x1
      error stop 1
    end if
    if (stats(2) < reject_injected) then
      print '(A,I0,A,I0)', 'recovery case expected nrej >= ', reject_injected, ', got ', stats(2)
      error stop 1
    end if
  end subroutine

  subroutine check_reject_limit_path()
    x0 = 0.0_wp
    x1 = 1.0_wp
    x = x0
    y = [1.0_wp]
    w = [1.0_wp]
    h0 = 0.25_wp
    eps = 1.0e-6_wp

    always_reject = .true.
    reject_remaining = 0
    call stiff3(n,fun,x,y,x1,jac,h0,eps,w,idid,stats=stats)

    if (idid /= -4) then
      print '(A,I0)', 'limit case expected idid=-4, got ', idid
      error stop 1
    end if
    if (stats(2) < 11) then
      print '(A,I0)', 'limit case expected at least 11 rejections, got ', stats(2)
      error stop 1
    end if
    if (stats(1) /= 0) then
      print '(A,I0)', 'limit case expected zero accepted steps, got ', stats(1)
      error stop 1
    end if
    if (abs(x - x0) > 2.0_wp*spacing(x0)) then
      print '(A,ES12.4,A,ES12.4)', 'limit case expected x near x0, got ', x, ' expected ', x0
      error stop 1
    end if
  end subroutine

  subroutine fun(n,y,f,ires)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    integer, intent(inout) :: ires

    f(1) = -y(1)
    if (always_reject) then
      ires = -1
    else if (reject_remaining > 0) then
      ires = -1
      reject_remaining = reject_remaining - 1
    end if
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -1.0_wp
  end subroutine

end program
