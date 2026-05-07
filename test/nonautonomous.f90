program nonautonomous

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 1, nprint = 100
  real(wp), parameter :: atol = 5.0e-6_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps, y_exact

  y = [0.0_wp]
  w = [1.0_wp]
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 1.0e-3_wp
  eps = 1.0e-8_wp

  call stiff3(n,fun,jac,dfdx,out,nprint,x0,x1,h0,eps,w,y)

  y_exact = 0.5_wp*(sin(x1) - cos(x1) + exp(-x1))
  if (abs(y(1)-y_exact) > atol) then
    error stop 1
  end if

contains

  subroutine fun(n,t,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -y(1) + sin(t)
  end subroutine

  subroutine jac(n,t,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -1.0_wp
  end subroutine

  subroutine dfdx(n,t,y,fx)
    integer, intent(in) :: n
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: fx(n)
    fx(1) = cos(t)
  end subroutine

  subroutine out(t,y,ih,qa)
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
  end subroutine

end program
