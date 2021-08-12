module robertson_rate

  use stiff3_solver, only: wp => stiff3_wp

  implicit none

contains

  !
  ! Robertson rate equations
  !
  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -4.e-2_wp*y(1)+1.e4_wp*y(2)*y(3)
    f(2) = -f(1)-3.e7_wp*y(2)**2
    f(3) = 3.e7_wp*y(2)**2
  end subroutine

  !
  ! Robertson rate equations - Jacobian
  !
  subroutine dfun(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -4.e-2_wp
    df(1,2) = 1.e4_wp*y(3)
    df(1,3) = 1.e4_wp*y(2)
    df(2,1) = 4.e-2_wp
    df(2,2) = -1.e4_wp*y(3)-6.e7_wp*y(2)
    df(2,3) = -1.e4_wp*y(2)
    df(3,1) = 0.0_wp
    df(3,2) = 6.e7_wp*y(2)
    df(3,3) = 0.0_wp
  end subroutine

end module

program main

  use stiff3_solver, only: wp => stiff3_wp, stiff3
  use robertson_rate, only: fun, dfun

  implicit none

  integer, parameter :: n = 3
  real(wp) :: y(n), w(n)
  real(wp) :: h0, eps, x0, x1

  integer, parameter :: nprint = 1

! initial value
  y = [1.0_wp, 0.0_wp, 0.0_wp]

! initial step size
  h0 = 0.0001_wp

! tolerance parameters
  eps = 1.e-5_wp
  w = 1
  w(2) = 1.e-4_wp

! time interval
  x0 = 0.0_wp
  x1 = 10.0_wp

  call output(x0,y,0,0.0_wp)
  call stiff3(n,fun,dfun,output,nprint,x0,x1,h0,eps,w,y)

contains

  subroutine output(x,y,iha,qa)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    real(wp) :: y2
    y2 = 1.e4_wp*y(2)
    print '(4(E19.12,2X),I4,2X,E19.12)', x, y(1), y2, y(3), iha, qa
  end subroutine

end program