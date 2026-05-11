program three_equation_system

  use stiff3_solver, only: wp => stiff3_wp, stiff3
  implicit none

  integer, parameter :: n = 3
  real(wp), parameter :: k1 = 0.013_wp
  real(wp), parameter :: k2 = 1000.0_wp
  real(wp), parameter :: k3 = 2500.0_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  integer :: irtrn, nsteps

! initial value
  y = [1.0_wp, 1.0_wp, 0.0_wp]
! initial step size
  h0 = 2.9e-4_wp
! tolerance
  eps = 1.0e-4_wp
  w = 1
! time interval
  x0 = 0.0_wp
  x1 = 50.0_wp

  nsteps = 0
  irtrn = 0
  call out(0,x0,x0,y,0,0.0_wp,irtrn)
  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,solout=out)

  print '(A,I0)', 'Accepted steps required: ', nsteps

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)

    f(1) = -k1*y(1) - k2*y(2)*y(3)
    f(2) = -k3*y(2)*y(3)
    f(3) = -k1*y(1) - k2*y(1)*y(3) - k3*y(2)*y(3)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -k1
    df(1,2) = -k2*y(3)
    df(1,3) = -k2*y(2)

    df(2,1) = 0.0_wp
    df(2,2) = -k3*y(3)
    df(2,3) = -k3*y(2)

    df(3,1) = -k1 - k2*y(3)
    df(3,2) = -k3*y(3)
    df(3,3) = -k2*y(1) - k3*y(2)
  end subroutine

  subroutine out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    nsteps = nr
  end subroutine

end program
