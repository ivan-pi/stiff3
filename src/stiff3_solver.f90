! Copyright 2021 Ivan Pribec
! SPDX-License-Identifier: Apache-2.0

!> Semi-implicit Runge-Kutta method of third order
!>
!> This is a modernized version of the code originally given in
!>
!>   Villadsen, J., & Michelsen, M. L. (1978). Solution of differential
!>   equation models by polynomial approximation. Prentice-Hall, Inc.,
!>   1978.
!>
module stiff3_solver

  use stiff3_linalg, only: lu, back

  implicit none
  private

  public :: stiff3, stiff3_wp
  public :: rhs_sub, jacobian_sub, output_sub

  !> Kind parameter for working precision of stiff3 reals
  integer, parameter :: stiff3_wp = kind(1.0d0)

  !> Working precision with short name for internal use
  integer, parameter :: wp = stiff3_wp

  abstract interface
    !> Function to evaluate the right-hand side of a system of ODEs.
    !> It is assumed the system of ODEs is autonomous, meaning that
    !> the independent variable x, does not appear explicitly.
    subroutine rhs_sub(n,y,f)
      import wp
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n)
      real(wp), intent(inout) :: f(n)
    end subroutine

    !> User supplied subprogram for evaluation of the Jacobian.
    subroutine jacobian_sub(n,y,df)
      import wp
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n)
      real(wp), intent(inout) :: df(n,n)
    end subroutine

    !> User supplied subprogram for output.
    subroutine output_sub(x,y,iha,qa)
      import wp
      real(wp), intent(in) :: x
        !! Current value of the independent variable
      real(wp), intent(in) :: y(:)
        !! Current value of the dependent variable vector
      integer, intent(in) :: iha
        !! Number of bisections (unsuccesful integrations) in
        !! the current step
      real(wp), intent(in) :: qa
        !! Step-length acceleration factor
    end subroutine

  end interface


contains

  ! TODO: Check if the the statement `h0 = h` should appear before
  !       or after exiting the routine.

  !> Semi-implicit Runge-Kutta integrator routine
  !
  subroutine stiff3(n,fun,dfun,out,nprint,x0,x1,h0,eps,w,y)
    integer, intent(in) :: n
      !! Number of equations to be integrated.
    procedure(rhs_sub) :: fun
      !! User supplied subprogram for function evaluation.
    procedure(jacobian_sub) :: dfun
      !! User supplied subprogram for evaluation of the Jacobian.
    procedure(output_sub) :: out
      !! User supplied subprogram for output.
    integer, intent(in) :: nprint
      !! Printing interval. For `nprint = k` the solution is only printed.
      !! at every kth step.
    real(wp), intent(in) :: x0, x1
      !! Limits of the independent variable between which the differential
      !! equation is solved.
    real(wp), intent(inout) :: h0
      !! Suggested initial half-step length. On exit `h0` contains suggested
      !! value of half-step length for continued integration beyond `x1`.
    real(wp), intent(in) :: eps, w(n)
      !! Tolerance parameters.
    real(wp), intent(inout) :: y(n)
      !! Vector of dependent variables at `x0`. On exit `y` is the vector of
      !! dependent variables at `x1`.

    real(wp), dimension(n) :: yk1, yk2, ya, yold, yold1, f, fold
      !! Workspace for solution vector and right-hand side
    real(wp), dimension(n,n) :: df, dfold
      !! Workspace for jacobian arrays
    integer :: ip(n)
      !! Workspace for the pivot array

    integer :: icon, iha, i, j, nout
    real(wp) :: x, h, e, es, q, qa

  ! icon = 0 except for last step which ends exactly at x1
    icon = 0

    nout = 0
    x = x0
    h = h0

    outer: do

    ! last step - or first step longer than interval

      if (x + 2.0_wp*h >= x1) then
        h = (x1 - x)/2.0_wp
        icon = 1
      end if

    ! other steps - limit to one quarter of remaining interval

      if ((icon == 0) .and. (x + 4.0_wp*h > x1)) then
        h = (x1 - x)/4.0_wp
      end if

    ! evaluate function and jacobian

      call fun(n,y,f)
      call dfun(n,y,df)

    ! keep values which are used in half-step integration

      do i = 1, n
        yold(i) = y(i)
        fold(i) = f(i)
        do j = 1, n
          dfold(i,j) = df(i,j)
        end do
      end do

    ! perform full integration step

      call sirk3(n,fun,ip,f,y,yk1,yk2,df,2*h)

      do i = 1, n
        ya(i) = y(i)
        y(i) = yold(i)
        f(i) = fold(i)
        do j = 1, n
          df(i,j) = dfold(i,j)
        end do
      end do

    ! full step finished, start half-step integration
    ! iha counts number of steplength bisections

      iha = -1
      inner: do
        iha = iha + 1

        call sirk3(n,fun,ip,f,y,yk1,yk2,df,h)
        call fun(n,y,f)
        call dfun(n,y,df)

        yold1 = y

        call sirk3(n,fun,ip,f,y,yk1,yk2,df,h)

      ! half step integration finished
      ! compute deviation and compare with tolerance

        e = 0.0_wp
        do i = 1, n
          es = w(i)*abs(ya(i)-y(i))/(1.0_wp+abs(y(i)))
          e = max(e,es)
        end do
        q = e/eps
        qa = (4.0_wp*q)**0.25_wp
        if (q <= 1.0_wp) then
          exit inner
        end if

      ! deviation too large- return to half-step with smaller h

        do i = 1, n
          ya(i) = yold1(i)
          y(i) = yold(i)
          f(i) = fold(i)
          do j = 1, n
            df(i,j) = dfold(i,j)
          end do
        end do

        h = h/2.0_wp
        icon = 0

      end do inner

    ! adjust y-vector

      do i = 1, n
        y(i) = y(i) + (y(i) - ya(i))/7.0_wp
      end do
      x = x + 2*h

    !  compute new stepsize

      qa = 1.0_wp/(qa+1.0e-10_wp)
      if (qa > 3.0_wp) qa = 3.0_wp
      h = qa*h

    ! perform output if appropriate

      nout = nout + 1
      if (mod(nout,nprint) == 0 .or. icon == 1) then
        call out(x,y,iha,qa)
      end if

    ! exit main loop

      if (icon == 1) then
        h0 = h
        return
      end if

    end do outer

  end subroutine


  !> Single-step semi-implicit integration
  !
  subroutine sirk3(n,fun,ipiv,f,y,yk1,yk2,df,h)
    integer, intent(in) :: n
      !! Size of the system of ODEs
    procedure(rhs_sub) :: fun
      !! Function to evaluate the right hand side
    integer, intent(inout) :: ipiv(n)
      !! Integer workspace used to store pivots in the LU factorization
    real(wp), intent(inout) :: f(n)
      !! On input, array of rhs values at beginning of step
    real(wp), intent(inout) :: y(n)
      !! On input contains the current approximation of the dependent variables.
      !! On output contains the approximation at the new time.
    real(wp), intent(inout) :: yk1(n),yk2(n)
      !! Real workspace arrays used in the implicit Runge-Kutta rule
    real(wp), intent(inout) :: df(n,n)
      !! On input contains the Jacobian values J,
      !! On output contains the factorized matrix (I - h a J) = LU
    real(wp), intent(in) :: h
      !! Step size of the independent variable

    integer :: i

    real(wp), parameter :: a  =  0.4358665215084589_wp
    real(wp), parameter :: r1 =  1.037609496131859_wp
    real(wp), parameter :: r2 =  0.8349304838526377_wp
    real(wp), parameter :: r3 = -0.6302020887244523_wp
    real(wp), parameter :: r4 = -0.2423378912600452_wp

    real(wp), parameter :: DF_TOL = 1.0e-12_wp
      !! Jacobian cutoff (elements smaller than DF_TOL are set to zero)

    !
    ! form matrix (I - h a J)
    !
    where (abs(df) > DF_TOL)
      df = -h*a*df
    elsewhere
      df = 0.0_wp
    end where

    do i = 1, n
      df(i,i) = df(i,i) + 1.0_wp
    end do

    !
    ! perform triangular decomposition and evaluate k1
    !
    call lu(df,ipiv)
    call back(df,f,ipiv)

    do i = 1, n
      yk1(i) = h*f(i)
      yk2(i) = y(i) + 0.75_wp * yk1(i)
    end do
    call fun(n,yk2,f)
    call back(df,f,ipiv)

    !
    ! evaluate k2
    !
    do i = 1, n
      yk2(i) = h*f(i)
      y(i) = y(i) + r1 * yk1(i) + r2 * yk2(i)
      yk2(i) = r3 * yk1(i) + r4 * yk2(i)
    end do

    !
    ! evaluate k3
    ! for convenience stored in yk2
    !
    call back(df,yk2,ipiv)
    do i = 1, n
      y(i) = y(i) + yk2(i)
    end do

  end subroutine

end module
