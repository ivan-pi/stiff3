! hosea_driver.f90 -- Test problems from Hosea & Shampine (1996)
!
! We test the performance of STIFF3 on the examples in
!
! Hosea, M.E. & Shampine, L.F. (1996). Analysis and implementation of TR-BDF2
! Applied Numerical Mathematics 20, pg. 21-37
!
! TODO: evaluate the solution and estimate errors in accordance
!       with Hosea and Shampine, section 7, page 33
!
program hosea_driver
  use stiff3_solver, only: wp => stiff3_wp, stiff3
  implicit none

  integer :: n, stats(6)
  real(wp) :: x, xend, eps, h0
  real(wp) :: y(3), w(3)

  eps = 0.005_wp
  w = 1

  ! Sample problem #1
  n = 3
  x = 0; xend = 12
  y(1:3) = [1,0,0]
  h0 = eps/10
  call stiff3(n,fun1,x,y,xend,jac1,h0,eps,w,&
    stats=stats)
  call print_stats("Problem #1", stats)

  ! Sample problem #2
  n = 3
  x = 0; xend = 50
  y(1:3) = [1,1,0]
  h0 = eps/10
  call stiff3(n,fun2,x,y,xend,jac2,h0,eps,w,&
    stats=stats)
  call print_stats("Problem #2 (D4)", stats)

  ! Sample problem #3
  n = 2
  x = 0; xend = 20
  y(1:2) = [real(wp) :: 0,0.25_wp]
  h0 = eps/10
  call stiff3(n,fun3,x,y,xend,jac3,h0,eps,w,&
    stats=stats)
  call print_stats("Problem #3 (van der Pol)", stats)

  ! Sample problem #4
  n = 3
  x = 0; xend = 0.4e7_wp
  y(1:3) = [1,0,0]
  h0 = eps/10
  call stiff3(n,fun4,x,y,xend,jac4,h0,eps,w,&
    stats=stats)
  call print_stats("Problem #4 (Robertson)", stats)

contains

  subroutine print_stats(problem,stats)
    character(len=*), intent(in) :: problem
    integer, intent(in) :: stats(6)
    write(*,'(A)') problem
    write(*,'("  nacc = ", I0)') stats(1)
    write(*,'("  nrej = ", I0)') stats(2)
    write(*,'("  nfev = ", I0)') stats(3)
    write(*,'("  njac = ", I0)') stats(4)
    write(*,'("  nlu  = ", I0)') stats(5)
    write(*,'("  nsol = ", I0)') stats(6)
  end subroutine

  !
  ! Sample problem #1
  !
  subroutine fun1(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    associate(t => y(3))
      f(1) = -500*y(1) + 500*cos(t) - sin(t)
      f(2) =   -1*y(2) + sin(t) + cos(t)
      f(3) = 1.0_wp
    end associate
  end subroutine
  subroutine jac1(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    associate(t=>y(3))
      df(1,1) = -500
      df(2,1) = 0
      df(3,1) = 0

      df(1,2) = 0
      df(2,2) = -1
      df(3,2) = 0

      df(1,3) = -500*sin(t) - cos(t)
      df(2,3) = cos(t) - sin(t)
      df(3,3) = 0
    end associate
  end subroutine

  !
  ! Sample problem #2, problem D4
  !
  subroutine fun2(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)

      f(1) = -0.013_wp*y(1) - 1000*y(1)*y(3)
      f(2) = -2500*y(2)*y(3)
      f(3) = -0.013_wp*y(1) - 1000*y(1)*y(3) - 2500*y(2)*y(3)

  end subroutine
  subroutine jac2(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -0.013_wp - 1000*y(3)
    df(2,1) = 0
    df(3,1) = -0.013_wp - 1000*y(3)

    df(1,2) = 0
    df(2,2) = -2500*y(3)
    df(3,2) = -2500*y(3)

    df(1,3) = -1000*y(1)
    df(2,3) = -2500*y(2)
    df(3,3) = -1000*y(1) - 2500*y(2)

  end subroutine

  !
  ! Sample problem #3, van der Pol
  !
  subroutine fun3(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    real(wp), parameter :: eps = 1

    f(1) = y(2)
    f(2) = eps*(1 - y(1)**2)*y(2) - y(1)

  end subroutine
  subroutine jac3(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    real(wp), parameter :: eps = 1

    df(1,1) = 0
    df(2,1) = 1

    df(1,2) = -2*eps*y(2)*y(1) - 1
    df(2,2) = eps*(1 - y(1)**2)

  end subroutine


  !
  ! Sample problem #4, scaled Robertson problem
  !
  subroutine fun4(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)

      f(1) = -0.04_wp*y(1) + 1.0e4_wp*y(2)*y(3)
      f(2) =  0.04_wp*y(1) - 1.0e4_wp*y(2)*y(3) - 3.0E7_wp*y(2)**2
      f(3) =                                      3.0E7_wp*y(2)**2

  end subroutine
  subroutine jac4(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -0.04_wp
    df(2,1) =  0.04_wp
    df(3,1) = 0

    df(1,2) =  1.0e4_wp*y(3)
    df(2,2) = -1.0e4_wp*y(3) - 6.0E7_wp*y(2)
    df(3,2) = 6.0E7_wp*y(2)

    df(1,3) =  1.0e4_wp*y(2)
    df(2,3) = -1.0e4_wp*y(2)
    df(3,3) = 0

  end subroutine

end program