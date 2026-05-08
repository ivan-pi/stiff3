!> Dense output example: van der Pol oscillator with uniform output grid.
!>
!> The solver takes variable-size steps, but contd3() is used inside the
!> solout callback to evaluate the solution at 10 equally-spaced sub-points
!> within each accepted step, producing smooth, uniformly-sampled output.
program dense_output

  use stiff3_solver, only: wp => stiff3_wp, stiff3, contd3

  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: mu = 3.0_wp  ! stiffness parameter
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  integer :: irtrn

  y   = [1.0_wp, 1.0_wp]
  h0  = 1.0e-3_wp
  eps = 1.0e-5_wp
  w   = 1.0_wp
  x0  = 0.0_wp
  x1  = 30.0_wp

  ! Print header
  write(*,'(A)') '# t   y1   y2'

  ! Print initial condition
  write(*,'(3(E18.12,2X))') x0, y(1), y(2)

  call stiff3(n, fun, jac, x0, x1, h0, eps, w, y, solout=out)

contains

  subroutine fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = y(2)
    f(2) = mu*(1.0_wp - y(1)**2)*y(2) - y(1)
  end subroutine

  subroutine jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = 0.0_wp
    df(1,2) = 1.0_wp
    df(2,1) = mu*(-2.0_wp*y(1))*y(2) - 1.0_wp
    df(2,2) = mu*(1.0_wp - y(1)**2)
  end subroutine

  subroutine out(nr, told, t, y, ih, qa, irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told, t, y(:), qa
    integer, intent(in) :: ih
    integer, intent(inout) :: irtrn

    integer, parameter :: nsub = 10
    integer :: k
    real(wp) :: targ

    ! Emit nsub uniformly-spaced sub-points, then the step endpoint
    do k = 1, nsub
      targ = told + (t - told) * real(k, wp) / real(nsub, wp)
      write(*,'(3(E18.12,2X))') targ, contd3(1, targ), contd3(2, targ)
    end do
  end subroutine

end program
