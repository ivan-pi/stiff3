! caroll_driver.f90 -- A set of of ODE test problems
!
! This driver implements the test problems used in
!
! Caroll, John (1989). A composite integration scheme for the numerical
! solution of systems of ordinary differential equations.
! Journal of Computational and Applied Mathematics 25, 1-13
!
!
program caroll_driver
  use stiff3_solver, only: wp => stiff3_wp, stiff3
  use caroll_test_problems

  implicit none

  integer :: n
  real(wp) :: x, xend, h0, eps, gerr
  real(wp) :: w(4), y(4)
  integer :: stats(6)

  w = 1
  eps = 0.001_wp

!
! P1
!
  n = 3
  x = 0; xend = 100
  y(1:3) = [0,0,0]
  h0 = eps/20
  gerr = 0

  call stiff3(n,p1_fun,x,y,xend,p1_jac,h0,eps,w,&
    stats=stats,solout=p1_out)
  call print_stats("P1",stats)
  print *, "Global error: ", gerr, gerr < eps

!
! P2
!

  n = 3
  x = 0; xend = 20
  y(1:3) = [1,1,0]
  h0 = eps/20
  gerr = 0

  call stiff3(n,p2_fun,x,y,xend,p2_jac,h0,eps,w,&
    stats=stats,solout=p2_out)
  call print_stats("P2",stats)
  print *, "Global error: ", gerr, gerr < eps

!
! P3
!

  n = 3
  x = 0; xend = 25
  y(1:3) = [25498.0_wp/1500.0_wp,-16499.0_wp/1500.0_wp,x]
  h0 = eps/20
  gerr = 0

  call stiff3(n,p3_fun,x,y,xend,p3_jac,h0,eps,w,&
    stats=stats,solout=p3_out)
  call print_stats("P3",stats)
  print *, "Global error: ", gerr, gerr < eps

!
! P4
!
  n = 3
  x = 0; xend = 20
  y(1:3) = [real(wp) :: 2,1,0]
  h0 = eps/20
  gerr = 0

  call stiff3(n,p4_fun,x,y,xend,p4_jac,h0,eps,w,&
    stats=stats,solout=p4_out)
  call print_stats("P4",stats)
  print *, "Global error: ", gerr, gerr < eps

!
! P5
!
  n = 4
  x = 0; xend = 20
  y(1:4) = -1.0_wp
  h0 = eps/20
  gerr = 0

  call stiff3(n,p5_fun,x,y,xend,p5_jac,h0,eps,w,&
    stats=stats,solout=p5_out)
  call print_stats("P5",stats)
  print *, "Global error: ", gerr, gerr < eps

!
! P6 (Robertson)
!
  n = 3
  x = 0; xend = 20
  y(1:3) = [1,0,0]
  h0 = eps/20
  call stiff3(n,p6_fun,x,y,xend,p6_jac,h0,eps,w,&
    stats=stats)
  call print_stats("P6",stats)

!
! P7
!
  n = 3
  x = 0; xend = 400
  y(1:3) = [0,0,0]
  h0 = eps/20
  call stiff3(n,p7_fun,x,y,xend,p7_jac,h0,eps,w,&
    stats=stats)
  call print_stats("P7",stats)

!
! P8
!
  n = 2
  x = 0; xend = 20
  y(1:2) = [1,0]
  h0 = eps/20
  call stiff3(n,p8_fun,x,y,xend,p8_jac,h0,eps,w,&
    stats=stats)
  call print_stats("P8",stats)

!
! P9
!
  n = 2
  x = 0; xend = 5
  y(1:2) = [10000.0_wp/10004.0_wp,1.0_wp]
  h0 = eps/20
  call stiff3(n,p9_fun,x,y,xend,p9_jac,h0,eps,w,&
    stats=stats)
  call print_stats("P9",stats)

!
! P10
!
  n = 2
  x = 0; xend = 20
  y(1:2) = [5,5]
  h0 = eps/20
  call stiff3(n,p10_fun,x,y,xend,p10_jac,h0,eps,w,&
    stats=stats)
  call print_stats("P10",stats)


! ---------------------

  !
  ! P7, Convergence study with respect to step size
  ! See Table 18, page 12, in Caroll (1989)
  !
  print '(A)', "P7 Convergence study"
  block
    integer :: m
    real(wp) :: x, xend, h0, y(3)
    n = 3
    eps = 0.1_wp
    do m = 3, 7
      x = 0; xend = 400
      y(1:3) = [0,0,0]
      h0 = 1.0_wp/real(2**m,wp)
      call stiff3(n,p7_fun,x,y,xend,p7_jac,h0,eps,w,&
        stats=stats,hmax=(h0))
      print *, y(1:3), stats
    end do
  end block

contains

  subroutine print_stats(problem,stats)
    character(len=*), intent(in) :: problem
    integer, intent(in) :: stats(6)

    print '(A)', problem
    print '(A,3(I0,1X))', 'accepted rejected nfev: ', stats(1), stats(2), stats(3)
    print '(A,3(I0,1X))', 'njev nlu nsol: ', stats(4), stats(5), stats(6)

  end subroutine

  subroutine p1_out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    associate(sol => p1_sol(t))
      gerr = max(gerr,maxval(abs(y - sol),dim=1))
    end associate

  end subroutine

  subroutine p2_out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    associate(sol => p2_sol(t))
      gerr = max(gerr,maxval(abs(y - sol),dim=1))
    end associate

  end subroutine

  subroutine p3_out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    associate(sol => p3_sol(t))
      gerr = max(gerr,maxval(abs(y - sol),dim=1))
    end associate

  end subroutine

  subroutine p4_out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    associate(sol => p4_sol(t))
      print *, t, abs(y - sol)
      gerr = max(gerr,maxval(abs(y - sol),dim=1))
    end associate

  end subroutine

  subroutine p5_out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    associate(sol => p5_sol(t))
      gerr = max(gerr,maxval(abs(y - sol),dim=1))
    end associate

  end subroutine
end program