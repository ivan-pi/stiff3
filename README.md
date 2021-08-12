# stiff3

`stiff3` is a Fortran subprogram for solving stiff autonomous systems of ordinary differential equations (ODE's) using a semi-implicit Runge-Kutta method with three steps (SIRK3). The `stiff3` source code was originally published in the work:

> Villadsen, J., & Michelsen, M. L. (1978). *Solution of differential equation models by polynomial approximation*. Prentice-Hall, Inc.

This repository provides a refactored version with a simplified procedural interface.

## Usage

Basic use of the solver is demonstrated using the [Van der Pol oscillator](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator):

```fortran
program vanpol

  use stiff3_solver, only: wp => stiff3_wp, stiff3
  implicit none

  integer, parameter :: n = 2, nout = 1
  real(wp), parameter :: mu = 1000.0_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps

! initial value
  y = [1.0_wp, 1.0_wp]
! initial step size
  h0 = 0.001_wp
! tolerance
  eps = 1.0e-4_wp
  w = 1
! time interval
  x0 = 0.0_wp
  x1 = 3000.0_wp
! output initial condition
  call out(x0,y,0,0.0_wp)
! integrate system of ODEs
  call stiff3(n,fun,jac,out,nout,x0,x1,h0,eps,w,y)

contains

  subroutine fun(n,y,f)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    f(1) = y(2)
    f(2) = mu*(1.0_wp - y(1)**2)*y(2) - y(1)
  end subroutine

  subroutine jac(n,y,df)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: df(n,n)
    df(1,1) = 0.0_wp
    df(1,2) = 1.0_wp
    df(2,1) = mu*y(2)*(-2*y(1)) - 1.0_wp
    df(2,2) = mu*(1.0_wp - y(1)**2)
  end subroutine

  subroutine out(t,y,ih,qa)
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    write(*,'(3(E18.12,2X),I4,2X,G0)') t, y(1), y(2), ih, qa
  end subroutine

end program
```

## Requirements

Minimal requirements include:
* a recent Fortran compiler
* BLAS and LAPACK libraries

## Method

The semi-implicit Runge-Kutta method used by `stiff3` was first published in

> Caillaud, J. B., & Padmanabhan, L. (1971). An improved semi-implicit Runge-Kutta method for stiff systems. The Chemical Engineering Journal, 2(4), 227-232. https://doi.org/10.1016/0300-9467(71)85001-3

The adaptive stepsize selection strategy is described in Villadsen & Michelsen (1978), Section 8.2.3, pages 314 - 317.

## Contributing

We look forward to all types of contributions. If you would like to propose additional features or submit a bug report please open a new issue.

For students interested in CSE, here are some contribution ideas:
- Support for banded or sparse Jacobian matrices
- Use BLAS kernels for vector operations
- Continuous (dense) output of variables
- Extend `stiff3` to non-autonomous systems of ODE's
- Advanced stepsize control settings
- Improve the repository and code documentation
- Write a tutorial on how to use `stiff3`