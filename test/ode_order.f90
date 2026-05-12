program ode_order

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: ln2 = log(2.0_wp)
  real(wp), parameter :: min_order = 3.0_wp
  real(wp), parameter :: hs(3) = [0.2_wp, 0.1_wp, 0.05_wp]
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  real(wp) :: err(3), p12, p23
  integer :: i, idid

  w = 1.0_wp
  ! Use a very large tolerance so the order check focuses on the
  ! Runge-Kutta discretization error instead of tolerance-driven bisections.
  eps = 1.0e6_wp

  do i = 1, size(hs)
    y = [1.0_wp, 0.0_wp]
    x0 = 0.0_wp
    x1 = hs(i)
    h0 = hs(i)/2.0_wp

    call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,idid)
    if (idid /= 0) then
      print '(A,I0,A,I0)', 'expected successful solve for sample ', i, ', got idid=', idid
      error stop 1
    end if
    err(i) = max(abs(y(1)-cos(x1)),abs(y(2)+sin(x1)))
  end do

  p12 = log(err(1)/err(2))/ln2
  p23 = log(err(2)/err(3))/ln2

  if (min(p12,p23) < min_order) then
    print '(A,2(F8.4,1X),A,F8.4)', 'observed order too low: p12 p23 = ', p12, p23, ' min required = ', min_order
    error stop 1
  end if

contains

  subroutine fun(n,y,f, ires)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    integer, intent(inout) :: ires
    f(1) = y(2)
    f(2) = -y(1)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = 0.0_wp
    df(1,2) = 1.0_wp
    df(2,1) = -1.0_wp
    df(2,2) = 0.0_wp
  end subroutine

end program
