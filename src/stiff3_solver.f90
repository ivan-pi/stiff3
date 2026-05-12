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

  public :: stiff3, stiff3_auto, stiff3_work, stiff3_wp, stiff3_interp
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
    subroutine output_sub(nr,xold,x,y,iha,qa,irtrn)
      import wp
      integer, intent(in) :: nr
        !! Number of the current grid point, starting at 1 for the initial value
      real(wp), intent(in) :: xold
        !! Previous value of the independent variable
      real(wp), intent(in) :: x
        !! Current value of the independent variable
      real(wp), intent(in) :: y(:)
        !! Current value of the dependent variable vector
      integer, intent(in) :: iha
        !! Number of bisections (unsuccessful integrations) in
        !! the current step
      real(wp), intent(in) :: qa
        !! Step-length acceleration factor
      integer, intent(inout) :: irtrn
        !! If set < 0 by the callback, integration is interrupted
    end subroutine

  end interface

  interface stiff3
    module procedure stiff3_auto
    module procedure stiff3_work
  end interface

  interface stiff3_interp
    module procedure stiff3_interp_component
    module procedure stiff3_interp_vector
  end interface


contains

  ! TODO: Check if the the statement `h0 = h` should appear before
  !       or after exiting the routine.

  !> Semi-implicit Runge-Kutta integrator routine
  !
  subroutine stiff3_auto(n,fun,x,y,xend,jac,h0,eps,w,idid,solout,stats,hmax)
    integer, intent(in) :: n
      !! Number of equations to be integrated.
    procedure(rhs_sub) :: fun
      !! User supplied subprogram for function evaluation.
    procedure(jacobian_sub) :: jac
      !! User supplied subprogram for evaluation of the Jacobian.
    real(wp), intent(inout) :: x
      !! Independent variable at the start of integration. On exit it is the
      !! current integration point (`xend` on successful completion).
    real(wp), intent(in) :: xend
      !! Limits of the independent variable between which the differential
      !! equation is solved.
    real(wp), intent(inout) :: h0
      !! Suggested initial half-step length. On exit `h0` contains suggested
      !! value of half-step length for continued integration beyond `x1`,
      !! or the suggested step size at an interrupted callback return.
    real(wp), intent(in) :: eps, w(n)
      !! Tolerance parameters.
    integer, intent(out) :: idid
      !! Solver exit flag: `0` success, `-1` LU factorization failure,
      !! `-2` user interrupt from `solout`, `-3` step-size underflow.
    real(wp), intent(inout) :: y(n)
      !! Vector of dependent variables at `x`. On exit `y` is the vector of
      !! dependent variables at the returned `x`.
    procedure(output_sub), optional :: solout
      !! User supplied subprogram for output.
    integer, intent(out), optional :: stats(6)
      !! Statistics array with `[nacc, nrej, nfev, njev, nlu, nsol]`.
    real(wp), intent(in), optional :: hmax
      !! Maximum absolute half-step size. If absent or zero, defaults to
      !! `abs(xend - x)`.

    real(wp), dimension(n) :: yk1, yk2, ya, yold, yold1, f, fold
      !! Workspace for solution vector and right-hand side
    real(wp), dimension(n,n) :: df, dfold
      !! Workspace for jacobian arrays
    integer :: ip(n)
      !! Workspace for the pivot array

    call stiff3_core(n,fun,x,y,xend,jac,h0,eps,w,idid, &
                     yk1=yk1,yk2=yk2,ya=ya,yold=yold,yold1=yold1, &
                     f=f,fold=fold,df=df,dfold=dfold,ip=ip, &
                     solout=solout,stats=stats,hmax=hmax)

  end subroutine


  !> Semi-implicit Runge-Kutta integrator routine with explicit workspace
  !
  subroutine stiff3_work(n,fun,x,y,xend,jac,h0,eps,w,rwork,iwork,idid,solout,stats,hmax)
    integer, intent(in) :: n
      !! Number of equations to be integrated.
    procedure(rhs_sub) :: fun
      !! User supplied subprogram for function evaluation.
    procedure(jacobian_sub) :: jac
      !! User supplied subprogram for evaluation of the Jacobian.
    real(wp), intent(inout) :: x
      !! Independent variable at the start of integration. On exit it is the
      !! current integration point (`xend` on successful completion).
    real(wp), intent(in) :: xend
      !! Limits of the independent variable between which the differential
      !! equation is solved.
    real(wp), intent(inout) :: h0
      !! Suggested initial half-step length. On exit `h0` contains suggested
      !! value of half-step length for continued integration beyond `x1`,
      !! or the suggested step size at an interrupted callback return.
    real(wp), intent(in) :: eps, w(n)
      !! Tolerance parameters.
    real(wp), intent(inout) :: y(n)
      !! Vector of dependent variables at `x`. On exit `y` is the vector of
      !! dependent variables at the returned `x`.
    real(wp), intent(inout) :: rwork(n*(7 + 2*n))
      !! Real workspace of size `n*(7 + 2*n)`.
    integer, intent(inout) :: iwork(n)
      !! Integer workspace of size `n`.
    integer, intent(out) :: idid
      !! Solver exit flag: `0` success, `-1` LU factorization failure,
      !! `-2` user interrupt from `solout`, `-3` step-size underflow.
    procedure(output_sub), optional :: solout
      !! User supplied subprogram for output.
    integer, intent(out), optional :: stats(6)
      !! Statistics array with `[nacc, nrej, nfev, njev, nlu, nsol]`.
    real(wp), intent(in), optional :: hmax
      !! Maximum absolute half-step size. If absent or zero, defaults to
      !! `abs(xend - x)`.

    call stiff3_core(n,fun,x,y,xend,jac,h0,eps,w,idid, &
                     yk1   = rwork(1), &
                     yk2   = rwork(n+1), &
                     ya    = rwork(2*n+1), &
                     yold  = rwork(3*n+1), &
                     yold1 = rwork(4*n+1), &
                     f     = rwork(5*n+1), &
                     fold  = rwork(6*n+1), &
                     df    = rwork(7*n+1), &
                     dfold = rwork(7*n+n*n+1), &
                     ip    = iwork(1), &
                     solout=solout,stats=stats,hmax=hmax)

  end subroutine


  !> Semi-implicit Runge-Kutta integrator routine core implementation
  !
  subroutine stiff3_core(n,fun,x,y,xend,jac,h0,eps,w,idid, &
                         yk1,yk2,ya,yold,yold1,f,fold,df,dfold,ip, &
                         solout,stats,hmax)
    integer, intent(in) :: n
    procedure(rhs_sub) :: fun
    procedure(jacobian_sub) :: jac
    real(wp), intent(inout) :: x
      !! Independent variable at the start of integration. On exit it is the
      !! current integration point (`xend` on successful completion).
    real(wp), intent(in) :: xend
    real(wp), intent(inout) :: h0
    real(wp), intent(in) :: eps, w(n)
    integer, intent(out) :: idid
    real(wp), intent(inout) :: y(n)
    real(wp), intent(inout) :: yk1(n), yk2(n), ya(n), yold(n), yold1(n), f(n), fold(n)
    real(wp), intent(inout) :: df(n,n), dfold(n,n)
    integer, intent(inout) :: ip(n)
    procedure(output_sub), optional :: solout
    integer, intent(out), optional :: stats(6)
    real(wp), intent(in), optional :: hmax

    integer :: icon, iha, i, j, irtrn, lu_info
    integer :: nfev, njev, nlu, nacc, nrej, nsol
    real(wp) :: x_current, xold, h, e, es, q, qa, hmax_used
    logical :: have_f

  ! icon = 0 except for last step which ends exactly at x1
    icon = 0

    x_current = x
    if (present(hmax)) then
      if (hmax < 0.0_wp) error stop 'stiff3: hmax must be a non-negative real value'
      if (hmax == 0.0_wp) then
        hmax_used = abs(xend - x_current)
      else
        hmax_used = min(hmax,abs(xend - x_current))
      end if
    else
      hmax_used = abs(xend - x_current)
    end if
    h = min(h0,hmax_used)
    nfev = 0
    njev = 0
    nlu = 0
    nacc = 0
    nrej = 0
    nsol = 0
    have_f = .false.
    idid = 0

    if (present(solout)) then
      irtrn = 0
      call solout(1,x_current,x_current,y,0,0.0_wp,irtrn)
      if (irtrn < 0) then
        idid = -2
        h0 = h
        x = x_current
        if (present(stats)) stats = [nacc, nrej, nfev, njev, nlu, nsol]
        return
      end if
    end if

    outer: do

    ! last step - or first step longer than interval

      if ((x_current + 2.0_wp*h >= xend) .and. ((xend - x_current)/2.0_wp <= hmax_used)) then
        h = (xend - x_current)/2.0_wp
        icon = 1
      end if

    ! other steps - limit to one quarter of remaining interval

      if ((icon == 0) .and. (x_current + 4.0_wp*h > xend)) then
        h = (xend - x_current)/4.0_wp
      end if

      h = min(h,hmax_used)

    ! evaluate function and jacobian

      if (.not. have_f) then
        ! On the first accepted step there is no saved endpoint rhs yet.
        call fun(n,y,f)
        nfev = nfev + 1
      end if
      call jac(n,y,df)
      njev = njev + 1

    ! keep values which are used in half-step integration

      do i = 1, n
        yold(i) = y(i)
        fold(i) = f(i)
        do j = 1, n
          dfold(i,j) = df(i,j)
        end do
      end do

    ! perform full integration step

      call sirk3(n,fun,ip,f,y,yk1,yk2,df,2*h,nfev,nlu,nsol,lu_info)
      if (lu_info /= 0) then
        idid = -1
        h0 = h
        x = x_current
        if (present(stats)) stats = [nacc, nrej, nfev, njev, nlu, nsol]
        return
      end if

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

        call sirk3(n,fun,ip,f,y,yk1,yk2,df,h,nfev,nlu,nsol,lu_info)
        if (lu_info /= 0) then
          idid = -1
          y = yold
          f = fold
          df = dfold
          h0 = h
          x = x_current
          if (present(stats)) stats = [nacc, nrej, nfev, njev, nlu, nsol]
          return
        end if
        call fun(n,y,f)
        nfev = nfev + 1
        call jac(n,y,df)
        njev = njev + 1

        yold1 = y

        call sirk3(n,fun,ip,f,y,yk1,yk2,df,h,nfev,nlu,nsol,lu_info)
        if (lu_info /= 0) then
          idid = -1
          y = yold
          f = fold
          df = dfold
          h0 = h
          x = x_current
          if (present(stats)) stats = [nacc, nrej, nfev, njev, nlu, nsol]
          return
        end if

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
        nrej = nrej + 1
        if (abs(h) <= spacing(x_current)) then
          idid = -3
          y = yold
          f = fold
          df = dfold
          h0 = h
          x = x_current
          if (present(stats)) stats = [nacc, nrej, nfev, njev, nlu, nsol]
          return
        end if

      end do inner

    ! adjust y-vector

      do i = 1, n
        y(i) = y(i) + (y(i) - ya(i))/7.0_wp
      end do
      xold = x_current
      x_current = x_current + 2*h

    ! evaluate rhs at the accepted-step end so it is ready for the next step
    ! and available in the explicit-workspace path through rwork(5*n+1:6*n)
      call fun(n,y,f)
      nfev = nfev + 1
      have_f = .true.

    !  compute new stepsize

      qa = min(1.0_wp/(qa + eps),3.0_wp)
      h = min(qa*h,hmax_used)

    ! perform output if appropriate

      nacc = nacc + 1
      if (present(solout)) then
        irtrn = 0
        call solout(nacc+1,xold,x_current,y,iha,qa,irtrn)
        if (irtrn < 0) then
          idid = -2
          h0 = h
          x = x_current
          if (present(stats)) stats = [nacc, nrej, nfev, njev, nlu, nsol]
          return
        end if
      end if

    ! exit main loop

      if (icon == 1) then
        idid = 0
        h0 = h
        x = x_current
        if (present(stats)) stats = [nacc, nrej, nfev, njev, nlu, nsol]
        return
      end if

    end do outer

  end subroutine


  !> Interpolate one component of the accepted solution over `[xold, x]`.
  !!
  !! This routine is intended for use from within a `solout` callback after
  !! calling `stiff3_work`, which provides the required accepted-step data in
  !! `rwork`.
  subroutine stiff3_interp_component(xold,x,y,rwork,xeval,idx,yeval)
    real(wp), intent(in) :: xold, x
    real(wp), intent(in) :: y(:), rwork(:)
    real(wp), intent(in) :: xeval
    integer, intent(in) :: idx
    real(wp), intent(out) :: yeval

    integer :: n
    real(wp) :: h, s
    real(wp) :: a1, a2, b1, b2

    n = size(y)

    if (idx < 1 .or. idx > n) error stop 'stiff3_interp_component: idx out of range'

    h = x - xold
    if (h == 0.0_wp) then
      yeval = y(idx)
      return
    end if

    s = max(0.0_wp,min(1.0_wp,(xeval - xold)/h))
    a1 = (1.0_wp + 2.0_wp*s)*(s - 1.0_wp)**2
    a2 = (3.0_wp - 2.0_wp*s)*s**2
    b1 = h*s*(s - 1.0_wp)**2
    b2 = h*(s - 1.0_wp)*s**2

    ! The rwork slices below must match the named associations in stiff3_work:
    ! yold -> rwork(3*n+1:4*n), fold -> rwork(6*n+1:7*n),
    ! f    -> rwork(5*n+1:6*n), holding the accepted-step end derivative.
    yeval = a1*rwork(3*n + idx) + a2*y(idx) + &
            b1*rwork(6*n + idx) + b2*rwork(5*n + idx)
  end subroutine


  !> Interpolate the accepted solution vector over `[xold, x]`.
  !!
  !! This routine is intended for use from within a `solout` callback after
  !! calling `stiff3_work`, which provides the required accepted-step data in
  !! `rwork`.
  subroutine stiff3_interp_vector(xold,x,y,rwork,xeval,yeval)
    real(wp), intent(in) :: xold, x
    real(wp), intent(in) :: y(:), rwork(:)
    real(wp), intent(in) :: xeval
    real(wp), intent(out) :: yeval(:)

    integer :: n
    real(wp) :: h, s
    real(wp) :: a1, a2, b1, b2

    n = size(y)

    h = x - xold
    if (h == 0.0_wp) then
      yeval = y
      return
    end if

    s = max(0.0_wp,min(1.0_wp,(xeval - xold)/h))
    a1 = (1.0_wp + 2.0_wp*s)*(s - 1.0_wp)**2
    a2 = (3.0_wp - 2.0_wp*s)*s**2
    b1 = h*s*(s - 1.0_wp)**2
    b2 = h*(s - 1.0_wp)*s**2

    ! The rwork slices below must match the named associations in stiff3_work:
    ! yold -> rwork(3*n+1:4*n), fold -> rwork(6*n+1:7*n),
    ! f    -> rwork(5*n+1:6*n), holding the accepted-step end derivative.
    yeval = a1*rwork(3*n+1:4*n) + a2*y + &
            b1*rwork(6*n+1:7*n) + b2*rwork(5*n+1:6*n)
  end subroutine


  !> Single-step semi-implicit integration
  !
  subroutine sirk3(n,fun,ipiv,f,y,yk1,yk2,df,h,nfev,nlu,nsol,info)
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
    integer, intent(inout) :: nfev
      !! Number of right-hand side evaluations
    integer, intent(inout) :: nlu
      !! Number of LU decompositions
    integer, intent(inout) :: nsol
      !! Number of linear-system solves (back substitutions)
    integer, intent(out) :: info
      !! LU factorization status from LAPACK `*getrf`

    integer :: i

    real(wp), parameter :: a  =  0.4358665215084589_wp
    real(wp), parameter :: r1 =  1.037609496131859_wp
    real(wp), parameter :: r2 =  0.8349304838526377_wp
    real(wp), parameter :: r3 = -0.6302020887244523_wp
    real(wp), parameter :: r4 = -0.2423378912600452_wp

    !
    ! form matrix (I - h a J)
    !
    df = -h*a*df
    do i = 1, n
      df(i,i) = df(i,i) + 1.0_wp
    end do

    !
    ! perform triangular decomposition and evaluate k1
    !
    call lu(df,ipiv,info)
    if (info /= 0) return
    nlu = nlu + 1
    call back(df,f,ipiv)
    nsol = nsol + 1

    do i = 1, n
      yk1(i) = h*f(i)
      yk2(i) = y(i) + 0.75_wp * yk1(i)
    end do
    call fun(n,yk2,f)
    nfev = nfev + 1
    call back(df,f,ipiv)
    nsol = nsol + 1

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
    nsol = nsol + 1
    do i = 1, n
      y(i) = y(i) + yk2(i)
    end do

  end subroutine

end module
