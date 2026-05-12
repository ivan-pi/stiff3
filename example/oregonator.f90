program oregonator

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 3
  integer, parameter :: ncases = 6
  real(wp), parameter :: y0(n) = [1.0_wp, 2.0_wp, 3.0_wp]
  real(wp), parameter :: yref(n) = [ &
    1.000814870318523_wp, 1228.178521549917_wp, 132.0554942846706_wp ]
  real(wp), parameter :: x0 = 0.0_wp
  real(wp), parameter :: x1 = 360.0_wp
  real(wp), parameter :: h0_initial = 1.0e-3_wp
  real(wp), parameter :: verification_tol = 1.0e-8_wp
  real(wp), parameter :: eps_values(ncases) = [ &
    1.0e-4_wp, 1.0e-5_wp, 1.0e-6_wp, 1.0e-7_wp, 1.0e-8_wp, 1.0e-9_wp ]
  character(len=*), parameter :: csv_file = 'oregonator_work_precision.csv'
  character(len=*), parameter :: svg_file = 'oregonator_work_precision.svg'
  character(len=*), parameter :: solution_file = 'oregonator_solution.dat'

  integer :: i
  integer :: solution_unit
  integer :: verification_stats(6)
  integer :: stats(6,ncases)
  real(wp) :: y(n), w(n), errors(ncases), cpu_times(ncases), verification_error, verification_cpu
  logical :: solution_output_enabled

  w = 1.0_wp
  solution_output_enabled = .false.

  call integrate_case(1.0e-8_wp, y, verification_stats, verification_cpu, solution_file)
  verification_error = max_relative_error(y, yref)
  if (verification_error > verification_tol) then
    print '(A,ES12.4,A,ES12.4)', &
      'reference verification failed: max relative error=', verification_error, &
      ' exceeds ', verification_tol
    error stop 1
  end if

  print '(A,ES12.4)', 'verified reference solution with max relative error ', verification_error
  print '(A,3(1X,ES24.16))', 'y(360) =', y
  print '(A,1X,ES12.4,1X,A)', 'verification cpu time:', verification_cpu, 's'
  print '(A,6(I0,1X))', 'verification stats [nacc nrej nfev njev nlu nsol]: ', verification_stats

  do i = 1, ncases
    call integrate_case(eps_values(i), y, stats(:,i), cpu_times(i))
    errors(i) = max_relative_error(y, yref)
  end do

  call write_csv(csv_file, eps_values, errors, cpu_times, stats)
  call write_svg(svg_file, eps_values, cpu_times)

  print '(A)', 'work-precision data:'
  print '(A)', '  eps         error      cpu[s]      nfev   njev   nlu   nsol'
  do i = 1, ncases
    print '(3(1X,ES10.2),4(1X,I6))', eps_values(i), errors(i), cpu_times(i), &
      stats(3,i), stats(4,i), stats(5,i), stats(6,i)
  end do
  print '(A,1X,A)', 'wrote', csv_file
  print '(A,1X,A)', 'wrote', svg_file
  print '(A,1X,A)', 'wrote', solution_file

contains

  subroutine integrate_case(eps, y, stats, cpu_time_seconds, solution_filename)
    real(wp), intent(in) :: eps
    real(wp), intent(out) :: y(n)
    integer, intent(out) :: stats(6)
    real(wp), intent(out) :: cpu_time_seconds
    character(len=*), intent(in), optional :: solution_filename

    integer :: idid, ios
    real(wp), parameter :: min_cpu_time = 1.0e-12_wp
    real(wp) :: cpu_start, cpu_end, h0
    logical :: write_solution

    y = y0
    h0 = h0_initial
    write_solution = present(solution_filename)
    if (write_solution) then
      open(newunit=solution_unit, file=solution_filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) error stop 'failed to open solution output file'
      write(solution_unit,'(A)') '# t y1 y2 y3'
      solution_output_enabled = .true.
    end if

    call cpu_time(cpu_start)
    call stiff3(n, fun, x0, y, x1, jac, h0, eps, w, idid, solout=output_solution, stats=stats)
    if (idid /= 0) then
      print '(A,I0)', 'stiff3 failed with idid=', idid
      error stop 1
    end if
    call cpu_time(cpu_end)
    cpu_time_seconds = max(min_cpu_time, cpu_end - cpu_start)

    if (write_solution) then
      close(solution_unit)
      solution_output_enabled = .false.
    end if
  end subroutine

  subroutine output_solution(nr, xold, x, y, iha, qa, irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: xold
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: iha
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn

    if (solution_output_enabled) then
      write(solution_unit,'(ES24.16,1X,*(ES24.16,1X))') x, y
    end if
  end subroutine

  function max_relative_error(y, y_reference) result(err)
    real(wp), intent(in) :: y(n), y_reference(n)
    real(wp) :: err

    err = maxval(abs(y - y_reference)/max(1.0_wp, abs(y_reference)))
  end function

  subroutine write_csv(filename, eps, err, cpu, stats)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: eps(:), err(:), cpu(:)
    integer, intent(in) :: stats(:,:)

    integer :: i, ios, unit

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open CSV output file'

    write(unit,'(A)') 'eps,max_relative_error,cpu_time_seconds,nfev,njev,nlu,nsol,nacc,nrej'
    do i = 1, size(eps)
      write(unit,'(ES24.16,'','',ES24.16,'','',ES24.16,'','',I0,'','',I0,'','',I0,'','',I0,'','',I0,'','',I0)') &
        eps(i), err(i), cpu(i), stats(3,i), stats(4,i), stats(5,i), stats(6,i), stats(1,i), stats(2,i)
    end do

    close(unit)
  end subroutine

  subroutine write_svg(filename, eps, cpu)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: eps(:), cpu(:)

    integer, parameter :: width = 800
    integer, parameter :: height = 600
    real(wp), parameter :: left = 100.0_wp
    real(wp), parameter :: right = 30.0_wp
    real(wp), parameter :: top = 50.0_wp
    real(wp), parameter :: bottom = 80.0_wp

    character(len=32) :: label
    integer :: i, ios, unit
    real(wp) :: log_cpu(size(cpu)), log_eps(size(eps))
    real(wp) :: xmin, xmax, ymin, ymax, xpad, ypad
    real(wp) :: x1p, x2p, y1p, y2p, xtick, ytick, tick_value

    log_cpu = log10(max(cpu, tiny(1.0_wp)))
    log_eps = log10(eps)

    xmin = minval(log_eps)
    xmax = maxval(log_eps)
    ymin = minval(log_cpu)
    ymax = maxval(log_cpu)
    xpad = 0.05_wp*max(xmax - xmin, 1.0_wp)
    ypad = 0.08_wp*max(ymax - ymin, 1.0_wp)
    xmin = xmin - xpad
    xmax = xmax + xpad
    ymin = ymin - ypad
    ymax = ymax + ypad

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open SVG output file'

    write(unit,'(A)') '<?xml version="1.0" encoding="UTF-8"?>'
    write(unit,'(A,I0,A,I0,A)') '<svg xmlns="http://www.w3.org/2000/svg" width="', width, &
      '" height="', height, '" viewBox="0 0 800 600">'
    write(unit,'(A)') '<rect width="100%" height="100%" fill="white"/>'
    write(unit,'(A)') &
      '<text x="400" y="28" font-family="sans-serif" font-size="22" ' // &
      'text-anchor="middle">Oregonator work-precision diagram</text>'
    write(unit,'(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)') &
      '<rect x="', left, '" y="', top, '" width="', real(width,wp) - left - right, &
      '" height="', real(height,wp) - top - bottom, '" fill="none" stroke="black" stroke-width="1.5"/>'

    do i = 0, 4
      xtick = left + real(i,wp)*(real(width,wp) - left - right)/4.0_wp
      tick_value = xmin + real(i,wp)*(xmax - xmin)/4.0_wp
      write(label,'(ES10.2)') 10.0_wp**tick_value
      write(unit,'(A,F0.3,A,F0.3,A)') '<line x1="', xtick, '" y1="520" x2="', xtick, &
        '" y2="526" stroke="#666"/>'
      write(unit,'(A,F0.3,A)') '<text x="', xtick, '" y="548" font-family="sans-serif" font-size="12" text-anchor="middle">'
      write(unit,'(A)') trim(adjustl(label))//'</text>'
    end do

    do i = 0, 4
      ytick = real(height,wp) - bottom - real(i,wp)*(real(height,wp) - top - bottom)/4.0_wp
      tick_value = 10.0_wp**(ymin + real(i,wp)*(ymax - ymin)/4.0_wp)
      write(label,'(ES10.2)') tick_value
      write(unit,'(A,F0.3,A,F0.3,A)') '<line x1="94" y1="', ytick, '" x2="100" y2="', ytick, &
        '" stroke="#666"/>'
      write(unit,'(A,F0.3,A)') '<text x="84" y="', ytick + 4.0_wp, '" font-family="sans-serif" font-size="12" text-anchor="end">'
      write(unit,'(A)') trim(adjustl(label))//'</text>'
    end do

    write(unit,'(A)') &
      '<text x="435" y="585" font-family="sans-serif" font-size="16" ' // &
      'text-anchor="middle">error tolerance eps (log scale)</text>'
    write(unit,'(A)') &
      '<text x="24" y="285" font-family="sans-serif" font-size="16" ' // &
      'text-anchor="middle" transform="rotate(-90 24 285)">' // &
      'CPU time (s) (log scale)</text>'

    do i = 1, size(cpu) - 1
      x1p = xcoord(log_eps(i), xmin, xmax, width, left, right)
      y1p = ycoord(log_cpu(i), ymin, ymax, height, top, bottom)
      x2p = xcoord(log_eps(i + 1), xmin, xmax, width, left, right)
      y2p = ycoord(log_cpu(i + 1), ymin, ymax, height, top, bottom)
      write(unit,'(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)') '<line x1="', x1p, '" y1="', y1p, &
        '" x2="', x2p, '" y2="', y2p, '" stroke="#1f77b4" stroke-width="2"/>'
    end do

    do i = 1, size(cpu)
      x1p = xcoord(log_eps(i), xmin, xmax, width, left, right)
      y1p = ycoord(log_cpu(i), ymin, ymax, height, top, bottom)
      write(label,'(ES9.1)') eps(i)
      write(unit,'(A,F0.3,A,F0.3,A)') '<circle cx="', x1p, '" cy="', y1p, '" r="4.5" fill="#d62728"/>'
      write(unit,'(A,F0.3,A,F0.3,A)') '<text x="', x1p + 8.0_wp, '" y="', y1p - 8.0_wp, &
        '" font-family="sans-serif" font-size="12">'
      write(unit,'(A)') 'eps='//trim(adjustl(label))//'</text>'
    end do

    write(unit,'(A)') '</svg>'

    close(unit)

  end subroutine

  function xcoord(value, xmin, xmax, width, left, right) result(xplot)
    real(wp), intent(in) :: value, xmin, xmax, left, right
    integer, intent(in) :: width
    real(wp) :: xplot

    xplot = left + (value - xmin)*(real(width,wp) - left - right)/(xmax - xmin)
  end function

  function ycoord(value, ymin, ymax, height, top, bottom) result(yplot)
    real(wp), intent(in) :: value, ymin, ymax, top, bottom
    integer, intent(in) :: height
    real(wp) :: yplot

    yplot = real(height,wp) - bottom - &
      (value - ymin)*(real(height,wp) - top - bottom)/(ymax - ymin)
  end function

  subroutine fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)

    f(1) = 77.27_wp*(y(2) + y(1)*(1.0_wp - 8.375e-6_wp*y(1) - y(2)))
    f(2) = (y(3) - (1.0_wp + y(1))*y(2))/77.27_wp
    f(3) = 0.161_wp*(y(1) - y(3))
  end subroutine

  subroutine jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    df = 0.0_wp
    df(1,1) = 77.27_wp*(1.0_wp - 2.0_wp*8.375e-6_wp*y(1) - y(2))
    df(1,2) = 77.27_wp*(1.0_wp - y(1))
    df(2,1) = -y(2)/77.27_wp
    df(2,2) = -(1.0_wp + y(1))/77.27_wp
    df(2,3) = 1.0_wp/77.27_wp
    df(3,1) = 0.161_wp
    df(3,3) = -0.161_wp
  end subroutine

end program
