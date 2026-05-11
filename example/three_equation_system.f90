program three_equation_system

  use stiff3_solver, only: wp => stiff3_wp, stiff3
  implicit none

  integer, parameter :: n = 3
  real(wp), parameter :: k1 = 0.013_wp
  real(wp), parameter :: k2 = 1000.0_wp
  real(wp), parameter :: k3 = 2500.0_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  integer :: stats(6)

! initial value
  y = [1.0_wp, 1.0_wp, 0.0_wp]
! initial step size
  h0 = 2.9e-4_wp
! tolerance
  eps = 1.0e-4_wp
  w = 1.0_wp
! time interval
  x0 = 0.0_wp
  x1 = 50.0_wp

  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,stats=stats)
  print '(A,3(I0,1X))', 'nfev njev nlu: ', stats(1), stats(2), stats(3)
  print '(A,3(I0,1X))', 'accepted rejected nsol: ', stats(4), stats(5), stats(6)

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)

    f(1) = -k1*y(1) - k2*y(1)*y(3)
    f(2) = -k3*y(2)*y(3)
    f(3) = -k1*y(1) - k2*y(1)*y(3) - k3*y(2)*y(3)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -k1 - k2*y(3)
    df(1,2) = 0.0_wp
    df(1,3) = -k2*y(1)

    df(2,1) = 0.0_wp
    df(2,2) = -k3*y(3)
    df(2,3) = -k3*y(2)

    df(3,1) = -k1 - k2*y(3)
    df(3,2) = -k3*y(3)
    df(3,3) = -k2*y(1) - k3*y(2)
  end subroutine

end program
