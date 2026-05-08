# stiff3

`stiff3` is a Fortran library for solving stiff autonomous systems of ordinary differential equations (ODEs) using a Rosenbrock-type (semi-implicit Runge-Kutta) method of third order. The `stiff3` source code was originally published in the work:

> Villadsen, J., & Michelsen, M. L. (1978). *Solution of differential equation models by polynomial approximation*. Prentice-Hall, Inc.

This repository provides a refactored version with a simplified procedural interface.

## Features

- Third-order Rosenbrock-type (semi-implicit Runge-Kutta) integrator suitable for stiff ODE systems
- Adaptive stepsize control with error tolerance per component
- Optional maximum absolute half-step size (`hmax`) to cap step growth
- Optional runtime statistics output: rhs evaluations (`nfev`), Jacobian evaluations (`njev`), and LU decompositions (`nlu`)
- Optional explicit workspace interface via `stiff3(..., rwork, iwork, ...)` for caller-managed memory
- Dense output helpers `stiff3_interp_component` and `stiff3_interp_vector` for accepted steps in `solout`
- Requires an exact user-supplied Jacobian
- Depends on BLAS and LAPACK for linear algebra operations
- Supports two build systems: [CMake](https://cmake.org/) and [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm)

## Installation

To use this project you need to have

* a recent Fortran compiler (e.g. GFortran 9+)
* BLAS and LAPACK libraries (e.g. OpenBLAS or the reference implementation)

### Using CMake

Configure, build, and run the tests with:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

To install the library to a custom prefix:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build
cmake --install build
```

### Using fpm

Install the [Fortran Package Manager](https://github.com/fortran-lang/fpm) and then build the library from the project root with:

```sh
fpm build
```

To use `stiff3` as a dependency in your own fpm project, add it to your `fpm.toml`:

```toml
[dependencies]
stiff3.git = "https://github.com/ivan-pi/stiff3"
```

### Running the examples

Four examples called `robertson`, `vanpol`, `lorenz`, and `predator_prey` are provided. They can be run with the command

```sh
fpm run --example <name>
```

With CMake, the compiled executables are placed in the `build/example/` directory and can be run directly or via CTest.

## Usage

Basic use of the solver is demonstrated using the [Van der Pol oscillator](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator):

```fortran
program vanpol

  use stiff3_solver, only: wp => stiff3_wp, stiff3
  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: mu = 10.0_wp
  real(wp) :: y(n), w(n), x0, x1, h0, eps, hmax
  integer :: irtrn, stats(3)

! initial value
  y = [1.0_wp, 1.0_wp]
! initial step size
  h0 = 0.001_wp
! tolerance
  eps = 1.0e-4_wp
  w = 1
! optional maximum absolute half-step size (0 means default abs(x1-x0))
  hmax = 0.0_wp
! time interval
  x0 = 0.0_wp
  x1 = 100.0_wp
! output initial condition
  irtrn = 0
  call out(0,x0,x0,y,0,0.0_wp,irtrn)
! integrate system of ODEs
  call stiff3(n,fun,x0,y,x1,jac,h0,eps,w,solout=out,stats=stats,hmax=hmax)
  print '(A,3(I0,1X))', 'nfev njev nlu: ', stats

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

  subroutine out(nr,told,t,y,ih,qa,irtrn)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: ih
    real(wp), intent(in) :: qa
    integer, intent(inout) :: irtrn
    write(*,'(3(E18.12,2X),I4,2X,G0)') t, y(1), y(2), ih, qa
  end subroutine

end program
```

## Method

`stiff3` implements a **Rosenbrock-type** (also known as semi-implicit Runge-Kutta) method. Rosenbrock methods linearize the implicit equations of a standard implicit Runge-Kutta scheme, requiring only a single LU factorization of the Jacobian per step. This makes them efficient for stiff problems while avoiding the nonlinear iterations of fully implicit methods.

The specific three-stage semi-implicit Runge-Kutta method (SIRK3) used by `stiff3` was first published in:

> Caillaud, J. B., & Padmanabhan, L. (1971). An improved semi-implicit Runge-Kutta method for stiff systems. *The Chemical Engineering Journal*, 2(4), 227–232. https://doi.org/10.1016/0300-9467(71)85001-3

The adaptive stepsize selection strategy is described in Villadsen & Michelsen (1978), Section 8.2.3, pages 314–317.

The two tolerance parameters `eps` and `w` together control how accurately the solution is tracked:

- **`eps`** — a single scalar that sets the overall accuracy goal. A smaller value requests a more accurate solution at the cost of more integration steps. Typical values range from `1e-3` (low accuracy) to `1e-8` (high accuracy).
- **`w(i)`** — a positive weight for each solution component `i`. Setting all weights to `1.0` applies the same accuracy goal uniformly across all components. Increasing `w(i)` relative to the others tightens the tolerance on component `i`.

A reasonable first choice is `eps = 1.0e-4` with all `w(i) = 1.0`, which requests roughly four significant digits from every component.

For interoperability scenarios where the caller manages memory allocation (e.g. language bindings), `stiff3` also provides an overload with explicit work arrays: `rwork` of size `n*(8 + 2*n)` and `iwork` of size `n`.

When using this explicit-workspace overload together with `solout`, accepted-step dense output is available through:

- `stiff3_interp_component(xold, x, y, rwork, xeval, icomp, yout)`
- `stiff3_interp_vector(xold, x, y, rwork, xeval, yout)`

These routines evaluate a cubic Hermite interpolant over the current accepted step `[xold, x]`. They are only valid inside the active `solout` callback of `stiff3_work`, and `xeval` must lie within the current accepted step.

You can also optionally pass `hmax` to limit the absolute half-step size used by the adaptive controller. If `hmax` is omitted or set to `0`, the default is `abs(x1-x0)`. If provided and positive, the solver uses `min(hmax, abs(x1-x0))`. Negative values are rejected.

**How the tolerance is applied:** after each step the solver estimates the local error in each component and checks:

```
|error_i| / (1 + |y_i|) <= eps / w(i)   for every i
```

The denominator `1 + |y_i|` makes the check relative when the solution is large and absolute when it is near zero, similar to the mixed-tolerance convention used in many modern ODE solvers. In standard `atol`/`rtol` notation this is equivalent to `atol_i = rtol_i = eps / w(i)` for each component.

## Contributing

We look forward to all types of contributions. If you would like to propose additional features or submit a bug report please open a new issue.

For students interested in CSE, here are some contribution ideas:
- Support for banded or sparse Jacobian matrices
- Use BLAS kernels for vector operations
- Continuous (dense) output of variables
- Extend `stiff3` to non-autonomous systems of ODEs
- Advanced stepsize control settings
- Write a tutorial on how to use `stiff3`
