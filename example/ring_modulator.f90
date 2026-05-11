program ring_modulator

  use stiff3_solver, only: stiff3, wp => stiff3_wp

  implicit none

  integer, parameter :: n = 16
  integer, parameter :: nphys = 15
  integer, parameter :: ncases = 33
  real(wp), parameter :: x0 = 0.0_wp
  real(wp), parameter :: x1 = 1.0e-3_wp
  real(wp), parameter :: h0_scale_factor = 1.0e2_wp
  real(wp), parameter :: verification_tol = 1.0e-6_wp
  real(wp), parameter :: y0(n) = 0.0_wp
  real(wp), parameter :: yref(nphys) = [ &
    -0.2339057358486745e-1_wp, &
    -0.7367485485540825e-2_wp, &
     0.2582956709291169e0_wp, &
    -0.4064465721283450e0_wp, &
    -0.4039455665149794e0_wp, &
     0.2607966765422943e0_wp, &
     0.1106761861269975e0_wp, &
     0.2939904342435596e-6_wp, &
    -0.2840029933642329e-7_wp, &
     0.7267198267264553e-3_wp, &
     0.7929487196960840e-3_wp, &
    -0.7255283495698965e-3_wp, &
    -0.7941401968526521e-3_wp, &
     0.7088495416976114e-4_wp, &
     0.2390059075236570e-4_wp ]
  character(len=*), parameter :: output_file = 'ring_modulator_work_precision.dat'

  integer :: i
  integer :: verification_stats(6)
  integer :: stats(6,ncases)
  real(wp) :: y(n), w(n), eps_values(ncases), errors(ncases), cpu_times(ncases), verification_error

  w = 1.0_wp
  do i = 1, ncases
    eps_values(i) = 10.0_wp**(-4.0_wp + real(i - 1, wp)/4.0_wp)
  end do

  call integrate_case(1.0e-10_wp, y, verification_stats, verification_error)
  if (verification_error > verification_tol) then
    print '(A,ES12.4,A,ES12.4)', &
      'reference verification failed: max relative error=', verification_error, &
      ' exceeds ', verification_tol
    error stop 1
  end if

  print '(A,ES12.4)', 'verified reference solution with max relative error ', verification_error
  print '(A,15(1X,ES24.16))', 'y(1.0e-3) =', y(1:nphys)
  print '(A,6(I0,1X))', 'verification stats [nacc nrej nfev njev nlu nsol]: ', verification_stats

  do i = 1, ncases
    call integrate_case(eps_values(i), y, stats(:,i), errors(i), cpu_times(i))
  end do

  call write_gnuplot_data(output_file, eps_values, cpu_times, errors, stats)

  print '(A)', 'work-precision data (CPU time vs eps):'
  print '(A)', '  eps         cpu[s]      error      nfev   njev   nlu   nsol'
  do i = 1, ncases
    print '(1X,ES10.2,1X,ES12.4,1X,ES10.2,4(1X,I10))', eps_values(i), cpu_times(i), errors(i), &
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
    real(wp) :: h0, cpu_start, cpu_end

    y = y0
    h0 = h0_scale_factor*eps
    call cpu_time(cpu_start)
    call stiff3(n, fun, x0, y, x1, jac, h0, eps, w, stats=stats)
    call cpu_time(cpu_end)
    if (present(err)) err = max_relative_error(y(1:nphys), yref)
    if (present(cpu_time_seconds)) cpu_time_seconds = max(min_cpu_time, cpu_end - cpu_start)
  end subroutine

  function max_relative_error(y, y_reference) result(err)
    real(wp), intent(in) :: y(nphys), y_reference(nphys)
    real(wp) :: err

    err = maxval(abs(y - y_reference)/max(1.0_wp, abs(y_reference)))
  end function

  subroutine write_gnuplot_data(filename, eps, cpu, err, stats)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: eps(:), cpu(:), err(:)
    integer, intent(in) :: stats(:,:)

    integer :: i, ios, unit

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open ring modulator output file'

    write(unit,'(A)') '# eps cpu_time_seconds max_relative_error nfev njev nlu nsol nacc nrej'
    do i = 1, size(eps)
      write(unit,'(ES24.16,1X,ES24.16,1X,ES24.16,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
        eps(i), cpu(i), err(i), stats(3,i), stats(4,i), stats(5,i), stats(6,i), stats(1,i), stats(2,i)
    end do

    close(unit)
  end subroutine

  subroutine fun(n, y, f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)

    real(wp), parameter :: c = 1.6e-8_wp
    real(wp), parameter :: cs = 2.0e-12_wp
    real(wp), parameter :: cp = 1.0e-8_wp
    real(wp), parameter :: r = 2.5e4_wp
    real(wp), parameter :: rp = 5.0e1_wp
    real(wp), parameter :: lh = 4.45_wp
    real(wp), parameter :: ls1 = 2.0e-3_wp
    real(wp), parameter :: ls2 = 5.0e-4_wp
    real(wp), parameter :: ls3 = 5.0e-4_wp
    real(wp), parameter :: rg1 = 36.3_wp
    real(wp), parameter :: rg2 = 17.3_wp
    real(wp), parameter :: rg3 = 17.3_wp
    real(wp), parameter :: ri = 5.0e1_wp
    real(wp), parameter :: rc = 6.0e2_wp
    real(wp), parameter :: gamma = 40.67286402e-9_wp
    real(wp), parameter :: delta = 17.7493332_wp
    real(wp), parameter :: pi = 3.141592653589793238462643383_wp
    real(wp), parameter :: exp_limit = 300.0_wp

    real(wp) :: t, uin1, uin2, ud1, ud2, ud3, ud4
    real(wp) :: qud1, qud2, qud3, qud4

    t = y(16)
    uin1 = 0.5_wp*sin(2.0e3_wp*pi*t)
    uin2 = 2.0_wp*sin(2.0e4_wp*pi*t)
    ud1 = y(3) - y(5) - y(7) - uin2
    ud2 = -y(4) + y(6) - y(7) - uin2
    ud3 = y(4) + y(5) + y(7) + uin2
    ud4 = -y(3) - y(6) + y(7) + uin2

    qud1 = gamma*(exp(min(delta*ud1, exp_limit)) - 1.0_wp)
    qud2 = gamma*(exp(min(delta*ud2, exp_limit)) - 1.0_wp)
    qud3 = gamma*(exp(min(delta*ud3, exp_limit)) - 1.0_wp)
    qud4 = gamma*(exp(min(delta*ud4, exp_limit)) - 1.0_wp)

    f(1) = (y(8) - 0.5_wp*y(10) + 0.5_wp*y(11) + y(14) - y(1)/r)/c
    f(2) = (y(9) - 0.5_wp*y(12) + 0.5_wp*y(13) + y(15) - y(2)/r)/c
    f(3) = (y(10) - qud1 + qud4)/cs
    f(4) = (-y(11) + qud2 - qud3)/cs
    f(5) = (y(12) + qud1 - qud3)/cs
    f(6) = (-y(13) - qud2 + qud4)/cs
    f(7) = (-y(7)/rp + qud1 + qud2 - qud3 - qud4)/cp
    f(8) = -y(1)/lh
    f(9) = -y(2)/lh
    f(10) = (0.5_wp*y(1) - y(3) - rg2*y(10))/ls2
    f(11) = (-0.5_wp*y(1) + y(4) - rg3*y(11))/ls3
    f(12) = (0.5_wp*y(2) - y(5) - rg2*y(12))/ls2
    f(13) = (-0.5_wp*y(2) + y(6) - rg3*y(13))/ls3
    f(14) = (-y(1) + uin1 - (ri + rg1)*y(14))/ls1
    f(15) = (-y(2) - (rc + rg1)*y(15))/ls1
    f(16) = 1.0_wp
  end subroutine

  subroutine jac(n, y, df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)

    real(wp), parameter :: c = 1.6e-8_wp
    real(wp), parameter :: cs = 2.0e-12_wp
    real(wp), parameter :: cp = 1.0e-8_wp
    real(wp), parameter :: r = 2.5e4_wp
    real(wp), parameter :: rp = 5.0e1_wp
    real(wp), parameter :: lh = 4.45_wp
    real(wp), parameter :: ls1 = 2.0e-3_wp
    real(wp), parameter :: ls2 = 5.0e-4_wp
    real(wp), parameter :: ls3 = 5.0e-4_wp
    real(wp), parameter :: rg1 = 36.3_wp
    real(wp), parameter :: rg2 = 17.3_wp
    real(wp), parameter :: rg3 = 17.3_wp
    real(wp), parameter :: ri = 5.0e1_wp
    real(wp), parameter :: rc = 6.0e2_wp
    real(wp), parameter :: gamma = 40.67286402e-9_wp
    real(wp), parameter :: delta = 17.7493332_wp
    real(wp), parameter :: pi = 3.141592653589793238462643383_wp
    real(wp), parameter :: exp_limit = 300.0_wp

    real(wp) :: t, uin2, duin1, duin2, ud1, ud2, ud3, ud4
    real(wp) :: qpud1, qpud2, qpud3, qpud4

    df = 0.0_wp

    t = y(16)
    uin2 = 2.0_wp*sin(2.0e4_wp*pi*t)
    duin1 = 1.0e3_wp*pi*cos(2.0e3_wp*pi*t)
    duin2 = 4.0e4_wp*pi*cos(2.0e4_wp*pi*t)
    ud1 = y(3) - y(5) - y(7) - uin2
    ud2 = -y(4) + y(6) - y(7) - uin2
    ud3 = y(4) + y(5) + y(7) + uin2
    ud4 = -y(3) - y(6) + y(7) + uin2

    qpud1 = gamma*delta*exp(min(delta*ud1, exp_limit))
    qpud2 = gamma*delta*exp(min(delta*ud2, exp_limit))
    qpud3 = gamma*delta*exp(min(delta*ud3, exp_limit))
    qpud4 = gamma*delta*exp(min(delta*ud4, exp_limit))

    df(1,1) = -1.0_wp/(c*r)
    df(1,8) = 1.0_wp/c
    df(1,10) = -0.5_wp/c
    df(1,11) = 0.5_wp/c
    df(1,14) = 1.0_wp/c
    df(2,2) = -1.0_wp/(c*r)
    df(2,9) = 1.0_wp/c
    df(2,12) = -0.5_wp/c
    df(2,13) = 0.5_wp/c
    df(2,15) = 1.0_wp/c
    df(3,3) = (-qpud1 - qpud4)/cs
    df(3,5) = qpud1/cs
    df(3,6) = -qpud4/cs
    df(3,7) = (qpud1 + qpud4)/cs
    df(3,10) = 1.0_wp/cs
    df(3,16) = (qpud1 + qpud4)*duin2/cs
    df(4,4) = (-qpud2 - qpud3)/cs
    df(4,5) = -qpud3/cs
    df(4,6) = qpud2/cs
    df(4,7) = (-qpud2 - qpud3)/cs
    df(4,11) = -1.0_wp/cs
    df(4,16) = (-qpud2 - qpud3)*duin2/cs
    df(5,3) = qpud1/cs
    df(5,4) = -qpud3/cs
    df(5,5) = (-qpud1 - qpud3)/cs
    df(5,7) = (-qpud1 - qpud3)/cs
    df(5,12) = 1.0_wp/cs
    df(5,16) = (-qpud1 - qpud3)*duin2/cs
    df(6,3) = -qpud4/cs
    df(6,4) = qpud2/cs
    df(6,6) = (-qpud2 - qpud4)/cs
    df(6,7) = (qpud2 + qpud4)/cs
    df(6,13) = -1.0_wp/cs
    df(6,16) = (qpud2 + qpud4)*duin2/cs
    df(7,3) = (qpud1 + qpud4)/cp
    df(7,4) = (-qpud2 - qpud3)/cp
    df(7,5) = (-qpud1 - qpud3)/cp
    df(7,6) = (qpud2 + qpud4)/cp
    df(7,7) = (-qpud1 - qpud2 - qpud3 - qpud4 - 1.0_wp/rp)/cp
    df(7,16) = (-qpud1 - qpud2 - qpud3 - qpud4)*duin2/cp
    df(8,1) = -1.0_wp/lh
    df(9,2) = -1.0_wp/lh
    df(10,1) = 0.5_wp/ls2
    df(10,3) = -1.0_wp/ls2
    df(10,10) = -rg2/ls2
    df(11,1) = -0.5_wp/ls3
    df(11,4) = 1.0_wp/ls3
    df(11,11) = -rg3/ls3
    df(12,2) = 0.5_wp/ls2
    df(12,5) = -1.0_wp/ls2
    df(12,12) = -rg2/ls2
    df(13,2) = -0.5_wp/ls3
    df(13,6) = 1.0_wp/ls3
    df(13,13) = -rg3/ls3
    df(14,1) = -1.0_wp/ls1
    df(14,14) = -(ri + rg1)/ls1
    df(14,16) = duin1/ls1
    df(15,2) = -1.0_wp/ls1
    df(15,15) = -(rc + rg1)/ls1
  end subroutine

end program
