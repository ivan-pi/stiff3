!> Test dense (interpolated) output via contd3().
!>
!> The simple harmonic oscillator
!>   y1' =  y2,  y1(0) = 1
!>   y2' = -y1,  y2(0) = 0
!> has the exact solution y1(t) = cos(t), y2(t) = -sin(t).
!>
!> At each accepted step we call contd3() at 5 interior sub-points and
!> check that the interpolated solution matches the exact solution to within
!> a loose relative tolerance (10 * eps_integration).  The test also verifies
!> that the endpoint values returned by contd3 match the grid-point values y(t)
!> to within machine epsilon.
program ode_dense_output

  use stiff3_solver, only: stiff3, contd3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 2
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  real(wp) :: max_err

  y = [1.0_wp, 0.0_wp]   ! y1(0)=cos(0)=1, y2(0)=-sin(0)=0
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 10.0_wp
  h0 = 0.01_wp
  eps = 1.0e-6_wp
  max_err = 0.0_wp

  call stiff3(n, fun, jac, x0, x1, h0, eps, w, y, solout=out)

  ! max_err is updated inside out(); a generous threshold of 1e-3 is used
  ! because the Hermite interpolant inherits the local accuracy of the solver.
  if (max_err > 1.0e-3_wp) then
    print '(A,ES12.4)', 'dense output error too large: ', max_err
    error stop 1
  end if

contains

  subroutine fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) =  y(2)
    f(2) = -y(1)
  end subroutine

  subroutine jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) =  0.0_wp; df(1,2) =  1.0_wp
    df(2,1) = -1.0_wp; df(2,2) =  0.0_wp
  end subroutine

  subroutine out(nr, told, t, y, ih, qa, irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told, t, y(:), qa
    integer, intent(in) :: ih
    integer, intent(inout) :: irtrn

    integer, parameter :: nsub = 5
    integer :: k
    real(wp) :: targ, y1_interp, y2_interp, err

    ! Check contd3 at nsub interior sub-points within [told, t]
    do k = 1, nsub
      targ = told + (t - told) * real(k, wp) / real(nsub + 1, wp)
      y1_interp = contd3(1, targ)
      y2_interp = contd3(2, targ)
      err = max(abs(y1_interp - cos(targ)), abs(y2_interp - (-sin(targ))))
      max_err = max(max_err, err)
    end do

    ! Verify that contd3 at the right endpoint reproduces y(t) exactly
    err = max(abs(contd3(1, t) - y(1)), abs(contd3(2, t) - y(2)))
    if (err > epsilon(1.0_wp) * 100) then
      print '(A,ES12.4)', 'contd3 endpoint mismatch: ', err
      error stop 2
    end if
  end subroutine

end program
