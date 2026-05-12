program hires

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 8
  integer, parameter :: ncases = 6
  real(wp), parameter :: x0 = 0.0_wp
  real(wp), parameter :: x1 = 321.8122_wp
  real(wp), parameter :: h0_initial = 1.0e-4_wp
  real(wp), parameter :: verification_tol = 1.0e-8_wp
  real(wp), parameter :: eps_values(ncases) = [ &
    1.0e-4_wp, 1.0e-5_wp, 1.0e-6_wp, 1.0e-7_wp, 1.0e-8_wp, 1.0e-9_wp ]
  real(wp), parameter :: y0(n) = [ &
    1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 5.7e-3_wp ]
  real(wp), parameter :: yref(n) = [ &
    0.7371312573325668e-3_wp, &
    0.1442485726316185e-3_wp, &
    0.5888729740967575e-4_wp, &
    0.1175651343283149e-2_wp, &
    0.2386356198831331e-2_wp, &
    0.6238968252742796e-2_wp, &
    0.2849998395185769e-2_wp, &
    0.2850001604814231e-2_wp ]
  character(len=*), parameter :: output_file = 'hires_work_precision.dat'

  integer :: i
  integer :: verification_stats(6)
  integer :: stats(6,ncases)
  real(wp) :: y(n), w(n), errors(ncases), cpu_times(ncases), verification_error

  w = 1.0_wp

  call integrate_case(1.0e-9_wp, y, verification_stats, verification_error)
  if (verification_error > verification_tol) then
    print '(A,ES12.4,A,ES12.4)', &
      'reference verification failed: max relative error=', verification_error, &
      ' exceeds ', verification_tol
    error stop 1
  end if

  print '(A,ES12.4)', 'verified reference solution with max relative error ', verification_error
  print '(A,8(1X,ES24.16))', 'y(321.8122) =', y
  print '(A,6(I0,1X))', 'verification stats [nacc nrej nfev njev nlu nsol]: ', verification_stats

  do i = 1, ncases
    call integrate_case(eps_values(i), y, stats(:,i), errors(i), cpu_times(i))
  end do

  call write_gnuplot_data(output_file, eps_values, cpu_times, errors, stats)

  print '(A)', 'work-precision data (CPU time vs eps):'
  print '(A)', '  eps         cpu[s]      error      nfev   njev   nlu   nsol'
  do i = 1, ncases
    print '(1X,ES10.2,1X,ES12.4,1X,ES10.2,4(1X,I6))', eps_values(i), cpu_times(i), errors(i), &
      stats(3,i), stats(4,i), stats(5,i), stats(6,i)
  end do
  print '(A,1X,A)', 'wrote', output_file

contains

  subroutine integrate_case(eps, y, stats, err, cpu_time_seconds)
    real(wp), intent(in) :: eps
    real(wp), intent(out) :: y(n)
    integer, intent(out) :: stats(6)
    real(wp), intent(out), optional :: err, cpu_time_seconds

    real(wp), parameter :: min_cpu_time = 1.0e-12_wp
    real(wp) :: h0, cpu_start, cpu_end, x
    integer :: idid

    y = y0
    x = x0
    h0 = h0_initial
    call cpu_time(cpu_start)
    call stiff3(n, fun, x, y, x1, jac, h0, eps, w, idid, stats=stats)
    call cpu_time(cpu_end)
    if (idid /= 0) error stop 'stiff3 failed in hires example'
    if (present(err)) err = max_relative_error(y, yref)
    if (present(cpu_time_seconds)) cpu_time_seconds = max(min_cpu_time, cpu_end - cpu_start)
  end subroutine

  function max_relative_error(y, y_reference) result(err)
    real(wp), intent(in) :: y(n), y_reference(n)
    real(wp) :: err

    err = maxval(abs(y - y_reference)/max(1.0_wp, abs(y_reference)))
  end function

  subroutine write_gnuplot_data(filename, eps, cpu, err, stats)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: eps(:), cpu(:), err(:)
    integer, intent(in) :: stats(:,:)

    integer :: i, ios, unit

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open HIRES output file'

    write(unit,'(A)') '# eps cpu_time_seconds max_relative_error nfev njev nlu nsol nacc nrej'
    do i = 1, size(eps)
      write(unit,'(ES24.16,1X,ES24.16,1X,ES24.16,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
        eps(i), cpu(i), err(i), stats(3,i), stats(4,i), stats(5,i), stats(6,i), stats(1,i), stats(2,i)
    end do

    close(unit)
  end subroutine

  subroutine fun(n, y, f, ires)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)

    integer, intent(inout) :: ires
    f(1) = -1.71_wp*y(1) + 0.43_wp*y(2) + 8.32_wp*y(3) + 0.0007_wp
    f(2) = 1.71_wp*y(1) - 8.75_wp*y(2)
    f(3) = -10.03_wp*y(3) + 0.43_wp*y(4) + 0.035_wp*y(5)
    f(4) = 8.32_wp*y(2) + 1.71_wp*y(3) - 1.12_wp*y(4)
    f(5) = -1.745_wp*y(5) + 0.43_wp*(y(6) + y(7))
    f(6) = -280.0_wp*y(6)*y(8) + 0.69_wp*y(4) + 1.71_wp*y(5) - 0.43_wp*y(6) + 0.69_wp*y(7)
    f(7) = 280.0_wp*y(6)*y(8) - 1.81_wp*y(7)
    f(8) = -f(7)
  end subroutine

  subroutine jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df = 0.0_wp
    df(1,1) = -1.71_wp
    df(1,2) = 0.43_wp
    df(1,3) = 8.32_wp
    df(2,1) = 1.71_wp
    df(2,2) = -8.75_wp
    df(3,3) = -10.03_wp
    df(3,4) = 0.43_wp
    df(3,5) = 0.035_wp
    df(4,2) = 8.32_wp
    df(4,3) = 1.71_wp
    df(4,4) = -1.12_wp
    df(5,5) = -1.745_wp
    df(5,6) = 0.43_wp
    df(5,7) = 0.43_wp
    df(6,4) = 0.69_wp
    df(6,5) = 1.71_wp
    df(6,6) = -280.0_wp*y(8) - 0.43_wp
    df(6,7) = 0.69_wp
    df(6,8) = -280.0_wp*y(6)
    df(7,6) = 280.0_wp*y(8)
    df(7,7) = -1.81_wp
    df(7,8) = 280.0_wp*y(6)
    df(8,6) = -280.0_wp*y(8)
    df(8,7) = 1.81_wp
    df(8,8) = -280.0_wp*y(6)
  end subroutine

end program
