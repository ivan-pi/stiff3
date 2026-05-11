program lorenz

  use stiff3_solver, only: wp => stiff3_wp, stiff3

  implicit none

  integer, parameter :: n = 3
  real(wp), parameter :: sigma = 10.0_wp, rho = 28.0_wp, beta = 8.0_wp/3.0_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  integer :: irtrn, idid

! initial value
  y = [1.0_wp, 1.0_wp, 1.0_wp]
! initial step size
  h0 = 1.0e-3_wp
! tolerance
  eps = 1.0e-5_wp
  w = 1
! time interval
  x0 = 0.0_wp
  x1 = 40.0_wp
! output initial condition
  irtrn = 0
  call out(0,x0,x0,y,0,0.0_wp,irtrn)
! integrate system of ODEs
  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,idid,solout=out)

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = sigma*(y(2) - y(1))
    f(2) = y(1)*(rho - y(3)) - y(2)
    f(3) = y(1)*y(2) - beta*y(3)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -sigma
    df(1,2) = sigma
    df(1,3) = 0.0_wp
    df(2,1) = rho - y(3)
    df(2,2) = -1.0_wp
    df(2,3) = -y(1)
    df(3,1) = y(2)
    df(3,2) = y(1)
    df(3,3) = -beta
  end subroutine

  subroutine out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn
    write(*,'(4(E18.12,2X),I4,2X,G0)') t, y(1), y(2), y(3), ih, qa
  end subroutine

end program
