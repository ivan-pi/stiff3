module ode_stats_problem

  use stiff3_solver, only: wp => stiff3_wp

  implicit none

  real(wp), parameter :: decay_rate = 2.0_wp
  integer :: fun_calls = 0
  integer :: jac_calls = 0

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    fun_calls = fun_calls + 1
    f(1) = -decay_rate*y(1)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    jac_calls = jac_calls + 1
    df(1,1) = -decay_rate
  end subroutine

end module

program ode_stats

  use stiff3_solver, only: stiff3, wp => stiff3_wp
  use ode_stats_problem, only: fun, jac, fun_calls, jac_calls

  implicit none

  integer, parameter :: n = 1
  real(wp) :: y(n), w(n), x0, x1, h0, eps
  integer :: idid, stats(6)

  y = [1.0_wp]
  w = 1.0_wp
  x0 = 0.0_wp
  x1 = 1.0_wp
  h0 = 0.01_wp
  eps = 1.0e-10_wp
  fun_calls = 0
  jac_calls = 0

  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,idid,stats=stats)

  associate(nacc => stats(1), nrej => stats(2), nfev => stats(3), njev => stats(4), nlu => stats(5), nsol => stats(6))
    if (nfev /= fun_calls .or. njev /= jac_calls) then
      print '(A,6(I0,1X),A,2(I0,1X))', &
        'stats mismatch [nacc nrej nfev njev nlu nsol]: ', stats, &
        ' callback counts [fun jac]: ', fun_calls, jac_calls
      error stop 1
    end if

    if (nlu <= 0) then
      print '(A,6(I0,1X))', 'expected positive nlu in stats [nacc nrej nfev njev nlu nsol], got: ', stats
      error stop 1
    end if

    if (nacc <= 0 .or. nrej < 0) then
      print '(A,6(I0,1X))', 'expected nonnegative step counters [nacc nrej nfev njev nlu nsol], got: ', stats
      error stop 1
    end if

    if (nsol < nlu) then
      print '(A,6(I0,1X))', 'expected nsol>=nlu in stats [nacc nrej nfev njev nlu nsol], got: ', stats
      error stop 1
    end if
  end associate

end program
