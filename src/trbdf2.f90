
module trbdf2_solver
use stiff3_solver, only: wp => stiff3_wp, &
    rhs_sub, jacobian_sub, output_sub
implicit none

! TRBDF2 coefficents
real(wp), parameter :: theta = 0.5_wp, &
    gamma = (1.0_wp - 1.0_wp/sqrt(2.0_wp))/theta

! BDF2 step coefficients
real(wp), parameter :: a2 = 2*(1-gamma*theta)/(1-2*gamma*theta)
real(wp), parameter :: a1 = (1 - a2)/gamma
real(wp), parameter :: a0 = -a1-a2

! other settings
real(wp), parameter :: beta = 2

integer, parameter :: itmax = 5
integer, parameter :: kjac = 15

contains

real(wp) function r_err(n,h,fold,fmid,fnew,err)
    integer, intent(in) :: n
    real(wp), intent(in) :: h, fold(n), fmid(n), fnew(n), err(n)
    real(wp) :: fac, tau, rsqr
    integer :: i

    ! coefficients in the error estimate
    real(wp), parameter :: t1 = 1/gamma
    real(wp), parameter :: t2 = 1/(gamma*(1-gamma))
    real(wp), parameter :: t3 = 1/(1-gamma)

    fac = (3*gamma**2*theta - 4*gamma*theta + 1)/(12*(1 - gamma*theta))

    ! Estimate error
    rsqr = 0
    do i = 1, n
        tau = abs(fac*h*2.0_wp*(t1*fold(i) - t2*fmid(i) + t3*fnew(i)))
        rsqr = rsqr + (tau/err(i))**2
    end do
    rsqr = rsqr / n
    r_err = sqrt(rsqr)

end function

subroutine trbdf2(n,fun,x,y,xend,jac,h,hmax,atol,rtol,stats,solout)
    integer, intent(in) :: n
    procedure(rhs_sub) :: fun
    procedure(jacobian_sub) :: jac
    real(wp), intent(inout) :: x, y(n)
    real(wp), intent(in) :: xend
    real(wp), intent(inout) :: h
    real(wp), value :: hmax
    real(wp), intent(in) :: rtol, atol(n)
    integer, intent(out) :: stats(6)
    procedure(output_sub), optional :: solout

    real(wp) :: B(n,n)
    real(wp), dimension(n) :: f,ymid,fmid,ynew,fnew,delta,err
    integer :: ipiv(n)

    integer :: nacc, nrej, nfev, njev, nlu, nsol
    logical :: new_jac, last_step, converged
    integer :: ijac, qstep, irtrn, iha

    real(wp) :: x_old, x_cur, r, q, hold

    ! PI Controller parameters
    real(wp), parameter :: kI = 0.233_wp
    real(wp), parameter :: kP = 0.133_wp
    real(wp), parameter :: safety = 0.9_wp
    real(wp), parameter :: facmax = 2.0_wp
    real(wp), parameter :: facmin = 0.2_wp
    real(wp) :: r_old, fac
    logical :: reject_last

    x_old = x
    x_cur = x

    if (hmax == 0.0_wp) hmax = abs(xend - x_cur)/4

    h = min(h, hmax)
    hold = h

    nacc = 0
    nrej = 0
    nfev = 0
    njev = 0
    nlu  = 0
    nsol = 0

    new_jac = .true.
    last_step = .false.

    q = 1.0_wp
    qstep = 1

    call fun(n,y,f); nfev = nfev + 1
    if (present(solout)) then
        irtrn = 0
        call solout(nacc+1,x_cur,x_cur,y,0,q,irtrn)
        if (irtrn < 0) then
            ! caller wants to terminate
            x = x_cur
            stats = [nacc,nrej,nfev,njev,nlu,nsol]
            write(*,'(A,I0)') "trbdf2: user requested termination with irtrn = ", irtrn
            return
        end if
    end if

    r_old = 1.0_wp      ! Assume error was exactly at tolerance initially
    reject_last = .false.

    controller: do

        !print *, "step = ", nacc + 1, ", h = ", h, ", qstep = ", qstep
        ! determine new h
        last_step = (x_cur + h) >= xend .and. (xend - x_cur) <= hmax

        ! clip in case of last step
        if (last_step) then
            h = xend - x_cur
            qstep = 0
        else if (h > hmax) then
            h = hmax
            qstep = 0
        end if

        iha = 0
        inner: do

            ! attempt step
            if (new_jac) ijac = 0
            call nonlin(n,fun,jac,new_jac,x_cur,h,y,f,rtol,atol,&
                B,ymid,fmid,ynew,fnew,delta,err,ipiv,&
                nfev,njev,nlu,nsol,&
                converged)

            if (converged) then
                err = rtol*abs(ynew) + atol
                r = r_err(n,h,f,fmid,fnew,err)
                !print *, "iha, r, h: ", iha, r, h
                if (r < 1.0_wp) then
                    ! converged
                    exit inner
                end if
            else
                ! Newton iteration diverged, assign a heavy penalty
                r = 4.0_wp
            end if

            ! --- STEP REJECTED ---
            ! Calculate a smooth reduction factor, but clamp it to facmin
            fac = safety * (1.0_wp / r)**(1.0_wp / 3.0_wp)
            fac = max(facmin, fac)
            reject_last = .true.
            h = h * fac

            ! step-halving
            !h = 0.5_wp * h
            new_jac = .true.
            qstep = 0

            if (abs(h) < spacing(x_cur)) error stop 99
            nrej = nrej + 1
            iha = iha + 1

        end do inner

        y = ynew
        f = fnew
        x_old = x_cur
        x_cur = x_cur + h
        hold = h
        ijac = ijac + 1

if (.true.) then
        ! Calculate PI inflation factor
        if (reject_last) then
            ! If we just recovered from a rejection, do NOT inflate.
            ! Just let the step size stabilize.
            fac = 1.0_wp
        else
            ! Standard PI control step
            fac = safety * (1.0_wp / r)**kI * (r_old / r)**kP
            ! Bound the growth to prevent Newton solver shock
            fac = min(facmax, fac)
        end if

        ! Update memory for the next step
        ! (Bound r_old away from absolute zero to prevent NaN in next step)
        r_old = max(1.0e-4_wp, r)
        reject_last = .false.
        ! compute new step-size for the NEXT iteration
        h = h * fac
else
        if (r > 0.5_wp) then
            qstep = qstep + 1
        else
            ! TODO: three-step rule
            if (qstep >= 3) then
                h = (1/r) ** (1.0_wp/3.0_wp) * hold
                qstep = 0
            else
                qstep = qstep + 1
            end if
        end if
end if
        ! compute new step-size ratio
        q = h / hold

        new_jac = r > 0.85_wp .or. q >= beta .or. ijac > 15
        !new_jac = .true.

        nacc = nacc + 1
        call fun(n,y,f); nfev = nfev + 1
        if (present(solout)) then
            irtrn = 0
            call solout(nacc+1,x_old,x_cur,y,iha,q,irtrn)
            if (irtrn < 0) then
                ! caller wants to terminate
                x = x_cur
                stats = [nacc,nrej,nfev,njev,nlu,nsol]
                write(*,'(A,I0)') "trbdf2: user requested termination with irtrn = ", irtrn
                return
            end if
        end if

        if (last_step) then
            x = x_cur
            stats = [nacc, nrej, nfev, njev, nlu, nsol]
            return
        end if

        !if (nacc > 20) stop
    end do controller

    error stop "trbdf2: unexpected error"

end subroutine trbdf2

subroutine nonlin(n,fun,jac,new_jac,x,h,yold,fold,rtol,atol,&
    B,ymid,fmid,ynew,fnew,delta,err,ipiv,nfev,njev,nlu,nsol,converged)
    use stiff3_linalg, only: lapack_getrf, back

    integer, intent(in) :: n
    procedure(rhs_sub) :: fun
    procedure(jacobian_sub) :: jac
    logical, intent(inout) :: new_jac
    real(wp), intent(in) :: x, h
    real(wp), intent(in) :: yold(n), fold(n)
    real(wp), intent(in) :: rtol, atol(n)

    ! Workspace
    real(wp), intent(inout) :: B(n,n)
    real(wp), intent(inout) :: ymid(n), fmid(n)
    real(wp), intent(inout) :: ynew(n), fnew(n)
    real(wp), intent(inout) :: delta(n)
    real(wp), intent(inout) :: err(n)
    integer, intent(inout) :: ipiv(n)
    integer, intent(inout) :: nfev, njev, nlu, nsol

    logical, intent(out) :: converged

    real(wp) :: fac, sold, smid, snew
    integer :: it, info, i

    if (new_jac) then
        call jac(n,yold,B); njev = njev + 1
        B = -(h*gamma*theta)*B
        do i = 1, n
            B(i,i) = B(i,i) + 1.0_wp
        end do
        ! Factorize matrix
        call lapack_getrf(n,n,B,n,ipiv,info); nlu = nlu + 1
        if (info /= 0) error stop 1
    end if

    !
    ! Modified Newton method
    !
    ! Stage 1
    ymid = yold
    fmid = fold
    converged = .false.
    do it = 1, itmax
        delta = yold - ymid + gamma*h*((1-theta)*fold + theta*fmid)
        call back(B,delta,ipiv,info); nsol = nsol + 1
        if (info /= 0) error stop 2
        err = rtol * abs(ymid) + atol
        ymid = ymid + delta
        call fun(n,ymid,fmid); nfev = nfev + 1
        if (all(abs(delta) < 0.1_wp*err)) then
            converged = .true.
            exit
        end if
    end do
    if (.not. converged) return


    ! Stage 2
    ynew = ymid
    fnew = fmid
    converged = .false.
    do it = 1, itmax
        delta = (h*fnew - (a0*yold + a1*ymid + a2*ynew))/a2
        call back(B,delta,ipiv,info); nsol = nsol + 1
        if (info /= 0) error stop 3
        err = rtol * abs(ynew) + atol
        ynew = ynew + delta
        call fun(n,ynew,fnew); nfev = nfev + 1
        if (all(abs(delta) < 0.1_wp*err)) then
            converged = .true.
            exit
        end if
    end do

end subroutine nonlin

end module trbdf2_solver



