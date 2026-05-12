program ode_robertson_invariants

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 3
  real(wp), parameter :: sum_tol = 1.0e-10_wp
  real(wp), parameter :: min_tol = -1.0e-14_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  real(wp) :: mass
  integer :: idid

  y = [1.0_wp, 0.0_wp, 0.0_wp]
  w = [1.0_wp, 1.0e-4_wp, 1.0_wp]
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 1.0e-4_wp
  eps = 1.0e-7_wp

  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,idid)
  if (idid /= 0) then
    print '(A,I0)', 'expected successful solve, got idid=', idid
    error stop 1
  end if

  mass = sum(y)
  if (abs(mass - 1.0_wp) > sum_tol) then
    print '(A,ES12.4,A,ES12.4)', 'robertson mass drift=', abs(mass-1.0_wp), ' exceeds tolerance=', sum_tol
    error stop 1
  end if

  if (minval(y) < min_tol) then
    print '(A,ES12.4)', 'robertson solution has negative component=', minval(y)
    error stop 1
  end if

contains

  subroutine fun(n,y,f, ires)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    integer, intent(inout) :: ires
    f(1) = -4.0e-2_wp*y(1) + 1.0e4_wp*y(2)*y(3)
    f(2) = -f(1) - 3.0e7_wp*y(2)**2
    f(3) = 3.0e7_wp*y(2)**2
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = -4.0e-2_wp
    df(1,2) = 1.0e4_wp*y(3)
    df(1,3) = 1.0e4_wp*y(2)
    df(2,1) = 4.0e-2_wp
    df(2,2) = -1.0e4_wp*y(3) - 6.0e7_wp*y(2)
    df(2,3) = -1.0e4_wp*y(2)
    df(3,1) = 0.0_wp
    df(3,2) = 6.0e7_wp*y(2)
    df(3,3) = 0.0_wp
  end subroutine

end program
