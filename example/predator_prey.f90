program predator_prey

  use stiff3_solver, only: wp => stiff3_wp, stiff3

  implicit none

  integer, parameter :: n = 2, nout = 10
  real(wp), parameter :: alpha = 1.1_wp, beta = 0.4_wp
  real(wp), parameter :: delta = 0.1_wp, gamma = 0.4_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps

! initial value
  y = [10.0_wp, 5.0_wp]
! initial step size
  h0 = 1.0e-3_wp
! tolerance
  eps = 1.0e-5_wp
  w = 1
! time interval
  x0 = 0.0_wp
  x1 = 40.0_wp
! output initial condition
  call out(x0,y,0,0.0_wp)
! integrate system of ODEs
  call stiff3(n,fun,jac,out,nout,x0,x1,h0,eps,w,y)

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = alpha*y(1) - beta*y(1)*y(2)
    f(2) = delta*y(1)*y(2) - gamma*y(2)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = alpha - beta*y(2)
    df(1,2) = -beta*y(1)
    df(2,1) = delta*y(2)
    df(2,2) = delta*y(1) - gamma
  end subroutine

  subroutine out(t,y,ih,qa)
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    write(*,'(3(E18.12,2X),I4,2X,G0)') t, y(1), y(2), ih, qa
  end subroutine

end program
