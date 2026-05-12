program pendant_drop

  use stiff3_solver, only: wp => stiff3_wp, stiff3

  implicit none

  integer, parameter :: n = 3
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  real(wp) :: pL, drho, ds
  integer :: stats(6), idid

! dimensionless pressure offset and density difference
  pL = 3.3888_wp
  drho = 0.8414_wp

! start from a short-distance asymptotic expansion to avoid y(2)=0 singularity
  ds = 1.0e-5_wp
  x0 = ds
  x1 = 2.0_wp

  y(1) = ds
  y(2) = cos(y(1))*ds
  y(3) = sin(y(1))*ds

! initial half-step length
  h0 = ds
! tolerance parameters
  eps = 1.0e-6_wp
  w = 1.0_wp

! integrate Young-Laplace system
  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,idid,solout=out,stats=stats)
  if (idid /= 0) error stop 'stiff3 failed in pendant_drop example'
  print '(A,3(I0,1X))', 'accepted rejected nfev: ', stats(1), stats(2), stats(3)
  print '(A,3(I0,1X))', 'njev nlu nsol: ', stats(4), stats(5), stats(6)

contains

  subroutine fun(n,y,f, ires)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    integer, intent(inout) :: ires
    f(1) = pL - drho*y(3) - sin(y(1))/y(2)
    f(2) = cos(y(1))
    f(3) = sin(y(1))
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -cos(y(1))/y(2)
    df(1,2) = sin(y(1))/y(2)**2
    df(1,3) = -drho
    df(2,1) = -sin(y(1))
    df(2,2) = 0.0_wp
    df(2,3) = 0.0_wp
    df(3,1) = cos(y(1))
    df(3,2) = 0.0_wp
    df(3,3) = 0.0_wp
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
