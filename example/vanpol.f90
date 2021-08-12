module vanpol

  use stiff3_solver, only: wp => stiff3_wp

  implicit none

  real(wp), parameter :: K = 10

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = y(2)
    f(2) = K*(1.0_wp - y(1)**2)*y(2) - y(1)
  end subroutine

  subroutine dfun(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = 0.0_wp
    df(1,2) = 1.0_wp
    df(2,1) = K*y(2)*(-2*y(1)) - 1.0_wp
    df(2,2) = K*(1.0_wp - y(1)**2)
  end subroutine

end module

program main

  use stiff3_solver, only: wp => stiff3_wp, stiff3
  use vanpol, only: fun, dfun

  implicit none

  integer, parameter :: n = 2
  real(wp) :: y(n), w(n)
  real(wp) :: h0, eps, x0, x1

  integer, parameter :: nprint = 1

! initial value
  y = [1.0_wp, 1.0_wp]

! initial step size
  h0 = 0.001_wp

! tolerance parameters
  eps = 1.e-4_wp
  w = 1

! time interval
  x0 = 0.0_wp
  x1 = 100.0_wp

  call output(x0,y,0,0.0_wp)
  call stiff3(n,fun,dfun,output,nprint,x0,x1,h0,eps,w,y)

contains

  subroutine output(x,y,iha,qa)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa

    print '(3(E18.12,2X),I4,2X,E18.12)', x, y(1), y(2), iha, qa
  end subroutine

end program