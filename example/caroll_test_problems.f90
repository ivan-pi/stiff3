module caroll_test_problems
implicit none

integer, parameter, private :: wp = kind(1.0d0)

contains

  ! P1
  function p1_sol(t) result(y)
    real(wp), intent(in) :: t
    real(wp) :: y(3)
    real(wp), parameter :: gamma = 94.0_wp/99.0_wp, lambda = 10001.0_wp
    y(1) = gamma*exp(-t) + &
      (10.0_wp/99.0_wp*exp(-1000*t) - 9496*cos(t) + 9506*sin(t))/lambda
    y(2) = gamma*exp(-t) + &
      (-108.0_wp/99.0_wp*exp(-1000*t) - 9494*cos(t)+9306*sin(t))/lambda
    y(3) = t
  end function

  subroutine p1_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    associate(t => y(3))
      f(1) = -6*y(1) + 5*y(2) + 2*sin(t)
      f(2) = 94*y(1) - 95*y(2)
      f(3) = 1.0_wp
    end associate
  end subroutine

  subroutine p1_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    associate(t=>y(3))
      df(1,1) = -6
      df(2,1) = 94
      df(3,1) = 0

      df(1,2) = 5
      df(2,2) = -95
      df(3,2) = 0

      df(1,3) = 2*cos(t)
      df(2,3) = 0
      df(3,3) = 0
    end associate
  end subroutine


  ! P2
  function p2_sol(t) result(y)
    real(wp), intent(in) :: t
    real(wp) :: y(3)
    y(1:2) = exp(-t)
    y(3) = t
  end function

  subroutine p2_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    real(wp), parameter :: alpha = 1, beta = 15
    associate(t => y(3))
      f(1) = -alpha*y(1) - beta*y(2) + (alpha + beta - 1) * exp(-t)
      f(2) = beta*y(1) - alpha*y(2) + (alpha - beta - 1) * exp(-t)
      f(3) = 1.0_wp
    end associate
  end subroutine

  subroutine p2_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    real(wp), parameter :: alpha = 1, beta = 15
    associate(t=>y(3))
      df(1,1) = -alpha
      df(2,1) = beta
      df(3,1) = 0

      df(1,2) = -beta
      df(2,2) = -alpha
      df(3,2) = 0

      df(1,3) = -(alpha+beta-1)*exp(-t)
      df(2,3) = -(alpha-beta-1)*exp(-t)
      df(3,3) = 0
    end associate
  end subroutine

  ! P3
  function p3_sol(t) result(y)
    real(wp), intent(in) :: t
    real(wp) :: y(3)
    y(1) = -2*exp(-t) + 7*exp(-1500*t) + (17998-14991*t)/1500
    y(2) = 1.5_wp*exp(-t) -3.5_wp*exp(-1500*t) - (13499-11245.5_wp*t)/1500
    y(3) = t
  end function
  subroutine p3_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    associate(t => y(3))
      f(1) = -4498*y(1) - 5996*y(2) + 0.006_wp - t
      f(2) = 2248.5_wp*y(1) + 2997*y(2) - 0.503_wp + 3*t
      f(3) = 1.0_wp
    end associate
  end subroutine
  subroutine p3_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    real(wp), parameter :: alpha = 1, beta = 15
    associate(t=>y(3))
      df(1,1) = -4498
      df(2,1) = 2248.5_wp
      df(3,1) = 0

      df(1,2) = -5996
      df(2,2) = 2997
      df(3,2) = 0

      df(1,3) = -1
      df(2,3) = 3
      df(3,3) = 0
    end associate
  end subroutine

  ! P4
  function p4_sol(t) result(y)
    real(wp), intent(in) :: t
    real(wp) :: y(3)
    y(1) = 2*F(t)
    y(2) = F(t)
    y(3) = t
  contains
    real(wp) function F(t)
      real(wp), intent(in) :: t
      real(wp), parameter :: b = 0.2_wp, mu = 1.0e-5_wp
      F = exp(-b*t) / (1 + mu*t)
    end function
  end function
  subroutine p4_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    real(wp), parameter :: b = 0.2_wp, g = 200, mu = 1.0e-5_wp

    associate(t => y(3))
      f(1) = -0.2_wp*((4*b+g)*y(1) + (2*b-2*g)*y(2)) - 2*mu/25*exp(-b*t)*(2*y(1) + y(2))**2
      f(2) = -0.2_wp*((2*b-2*g)*y(1) + (b+4*g)*y(2)) - mu/25*exp(-b*t)*(2*y(1) + y(2))**2
      f(3) = 1.0_wp
    end associate
  end subroutine

  subroutine p4_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    real(wp), parameter :: b = 0.2_wp, g = 200, mu = 1.0e-5_wp

    associate(t=>y(3))
      df(1,1) = -0.2_wp*(4*b+g) - 8*mu/25*exp(-b*t)*(2*y(1) + y(2))
      df(2,1) = -0.2_wp*(2*b-2*g) - 4*mu/25*exp(-b*t)*(2*y(1) + y(2))
      df(3,1) = 0

      df(1,2) = -0.2_wp*(2*b-2*g) - 4*mu/25*exp(-b*t)*(2*y(1) + y(2))
      df(2,2) = -0.2_wp*(b+4*g) - 2*mu/25*exp(-b*t)*(2*y(1) + y(2))
      df(3,2) = 0

      df(1,3) = b*2*mu/25*exp(-b*t)*(2*y(1) + y(2))**2
      df(2,3) = b*mu/25*exp(-b*t)*(2*y(1) + y(2))**2
      df(3,3) = 0
    end associate
  end subroutine

  ! P5
  function p5_sol(t) result(y)
    real(wp), intent(in) :: t
    real(wp) :: y(4)
    real(wp), parameter :: beta(4) = [1000.0_wp,800.0_wp,-10.0_wp,0.001_wp]
    y = beta / (1.0_wp - (1.0_wp + beta)*exp(beta*t))
  end function
  subroutine p5_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    real(wp), parameter :: beta(4) = [1000.0_wp,800.0_wp,-10.0_wp,0.001_wp]
    f(1:4) = -beta*y(1:4) - y(1:4)**2
  end subroutine
  subroutine p5_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    real(wp), parameter :: beta(4) = [1000.0_wp,800.0_wp,-10.0_wp,0.001_wp]
    integer :: i
    df(1:4,1:4) = 0
    do i = 1, 4
      df(i,i) = -beta(i) - 2*y(i)
    end do
  end subroutine

  ! P6
  subroutine p6_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -0.04_wp*y(1) + 1.0E4_wp*y(2)*y(3)
    f(2) =  0.04_wp*y(1) - 1.0E4_wp*y(2)*y(3) - 3.0E7_wp*y(2)**2
    f(3) =                                      3.0E7_wp*y(2)**2
  end subroutine
  subroutine p6_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -0.04_wp
    df(2,1) =  0.04_wp
    df(3,1) = 0

    df(1,2) =  1.0E4_wp*y(3)
    df(2,2) = -1.0E4_wp*y(3) - 6.0E7_wp*y(2)
    df(3,2) =                  6.0E7_wp*y(2)

    df(1,3) =  1.0E4_wp*y(2)
    df(2,3) = -1.0E4_wp*y(2)
    df(3,3) = 0

  end subroutine

  ! P7
  subroutine p7_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    associate(t => y(3))
      f(1) = 0.2_wp*(y(2) - y(1))
      f(2) = 10*y(1) - (60 - 0.125_wp*t)*y(2) + 0.125_wp*t
      f(3) = 1.0_wp
    end associate
  end subroutine
  subroutine p7_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    associate(t=>y(3))
      df(1,1) = -0.2_wp
      df(2,1) = 10
      df(3,1) = 0

      df(1,2) = 0.2_wp
      df(2,2) = 0.125_wp*t - 60
      df(3,2) = 0

      df(1,3) = 0
      df(2,3) = 0.125_wp
      df(3,3) = 0
    end associate
  end subroutine

  ! P8
  subroutine p8_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -y(2) + (1 - y(1)**2 - y(2)**2)
    f(2) =  y(1) + (1 - y(1)**2 - y(2)**2)
  end subroutine
  subroutine p8_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -2*y(1)
    df(2,1) = 1 - 2*y(1)

    df(1,2) = -1 - 2*y(2)
    df(2,2) = -2*y(2)

  end subroutine


  ! P9
  subroutine p9_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -10004*y(1) + 1.0E5_wp*y(2)**4
    f(2) = y(1) - y(2) - y(2)**4
  end subroutine

  subroutine p9_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -10004
    df(2,1) = 1

    df(1,2) = 4.0E5_wp*y(2)**3
    df(2,2) = -1 - 4*y(2)**3

  end subroutine


  ! P10
  subroutine p10_fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = -y(1)
    f(2) = y(1)**2 - 2*y(2)
  end subroutine

  subroutine p10_jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df(1,1) = -1
    df(2,1) = 2*y(1)

    df(1,2) = 0
    df(2,2) = -2

  end subroutine

end module