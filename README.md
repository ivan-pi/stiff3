# stiff3

**`stiff3` is a Fortran library for solving stiff, autonomous, first-order ODE systems using a 3rd-order Rosenbrock (semi-implicit Runge-Kutta) method.**

It solves initial value problems of the form:
```text
  y'(t) = f(y),   y(t0) = y0
```

This repository provides a heavily refactored version of the original implementation by Villadsen & Michelsen (1978).

> **Original source:** Villadsen, J., and Michelsen, M. L. *Solution of Differential Equation Models by Polynomial Approximation.* Prentice-Hall, Inc., 1978.

## Features

- Third-order Rosenbrock-type (semi-implicit Runge-Kutta) integrator suitable for stiff ODE systems
- Adaptive stepsize control with error tolerance per component
- Required integer exit flag (`idid`) reporting success and early-exit conditions
- Optional maximum absolute half-step size (`hmax`) to cap step growth
- Optional runtime statistics output: successful steps (`nacc`), rejected steps (`nrej`), rhs evaluations (`nfev`), Jacobian evaluations (`njev`), LU decompositions (`nlu`), and linear-system solves (`nsol`)
- Optional explicit workspace interface via `stiff3(..., rwork, iwork, ...)` for caller-managed memory
- Dense output helper `stiff3_interp` for accepted steps in `solout`
- Requires an exact user-supplied Jacobian
- Depends on BLAS and LAPACK for linear algebra operations
- Supports two build systems: [CMake](https://cmake.org/) and [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm)

## Releases

See all releases on GitHub: <https://github.com/ivan-pi/stiff3/releases>

| Major.Minor | Date       | Highlights |
|-------------|------------|------------|
| 0.2         | TBD        | Standardized solver call signature, optional runtime stats, optional explicit work arrays |
| 0.1         | 2023-02-13 | First public release |

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

Ten examples called `robertson`, `vanpol`, `lorenz`, `predator_prey`, `pendant_drop`, `three_equation_system`, `oregonator`, `e5_chemical`, `hires`, and `ring_modulator` are provided. They can be run with the command

```sh
fpm run --example <name>
```

With CMake, the compiled executables are placed in the `build/example/` directory and can be run directly or via CTest.

The `oregonator` example verifies a published reference solution at `t = 360` and writes `oregonator_work_precision.csv`, an `oregonator_work_precision.svg` plot of CPU time vs error tolerance `eps` (logarithmic y-axis), and `oregonator_solution.dat` (columns `t y1 y2 y3` for gnuplot) in the current working directory.

The `e5_chemical` example verifies the Datta/Bari E5 reference solution at `t = 1.0e13` and writes `e5_chemical_work_precision.csv` together with an `e5_chemical_work_precision.svg` plot of CPU time vs error tolerance `eps` (logarithmic y-axis).

The `hires` example verifies the Bari HIRES reference solution at `t = 321.8122` and writes `hires_work_precision.dat` (columns: `eps cpu_time_seconds max_relative_error nfev njev nlu nsol nacc nrej`) for plotting CPU time vs error tolerance `eps` with gnuplot.

The `ring_modulator` example verifies the Bari Ring Modulator reference solution at `t = 1.0e-3` and writes `ring_modulator_work_precision.dat` (columns: `eps cpu_time_seconds max_relative_error nfev njev nlu nsol nacc nrej`) for plotting CPU time vs error tolerance `eps` with gnuplot.

## Usage

Basic use of the solver is demonstrated using the [Van der Pol oscillator](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator):

```fortran
program vanpol

  use stiff3_solver, only: wp => stiff3_wp, stiff3
  implicit none

  integer, parameter :: n = 2
  real(wp), parameter :: mu = 10.0_wp
  real(wp) :: y(n), w(n), x0, x, x1, h0, eps, hmax
  integer :: idid, stats(6)

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
  x = x0
  x1 = 100.0_wp
! integrate system of ODEs
  call stiff3(n,fun,x,y,x1,jac,h0,eps,w,idid,solout=out,stats=stats,hmax=hmax)
  if (idid /= 0) error stop 'stiff3 failed'
! on successful exit, idid == 0 and x == x1
  print '(A)', 'Solver statistics'
  print '(A,I0)', '  nacc: ', stats(1)
  print '(A,I0)', '  nrej: ', stats(2)
  print '(A,I0)', '  nfev: ', stats(3)
  print '(A,I0)', '  njev: ', stats(4)
  print '(A,I0)', '  nlu: ', stats(5)
  print '(A,I0)', '  nsol: ', stats(6)

contains

  subroutine fun(n,y,f, ires)
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n)
    real(wp), intent(inout) :: f(n)
    integer, intent(inout) :: ires
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

For interoperability scenarios where the caller manages memory allocation (e.g. language bindings), `stiff3` also provides an overload with explicit work arrays: `rwork` of size `n*(7 + 2*n)` and `iwork` of size `n`. The required `idid` argument appears before the optional arguments in both overloads:

- `call stiff3(n, fun, x, y, x1, jac, h0, eps, w, idid, [solout=...], [stats=...], [hmax=...])`
- `call stiff3(n, fun, x, y, x1, jac, h0, eps, w, rwork, iwork, idid, [solout=...], [stats=...], [hmax=...])`

The integer exit flag `idid` reports solver status:

- `0` — successful completion at `x1`
- `-1` — LU factorization failed because the Jacobian system matrix was singular
- `-2` — integration interrupted by a user callback (`ires < -1` in `fun` or `irtrn < 0` in `solout`)
- `-3` — step-size underflow occurred during bisection (`abs(h) <= spacing(x_current)`)
- `-4` — too many physical rejections signaled by `fun` (`ires = -1`, limit exceeded)

The `fun` callback uses a DASSL-like `ires` contract:

- `ires = 0` — normal return
- `ires = -1` — reject the current trial step and retry with halved step size
- `ires < -1` — interrupt integration and return with `idid = -2`

When using this explicit-workspace overload together with `solout`, accepted-step dense output is available through the generic interface `stiff3_interp`:

- `stiff3_interp(xold, x, y, rwork, xeval, idx, yeval)`
- `stiff3_interp(xold, x, y, rwork, xeval, yeval)`

These routines evaluate a cubic Hermite interpolant over the current accepted step `[xold, x]`. They are only valid inside the active `solout` callback of `stiff3_work`, and `xeval` must lie within the current accepted step.

You can also optionally pass `hmax` to limit the absolute half-step size used by the adaptive controller. If `hmax` is omitted or set to `0`, the default is `abs(x1-x)` using the initial value of `x` at the start of integration. If provided and positive, the solver uses `min(hmax, abs(x1-x))`. Negative values are rejected.

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
