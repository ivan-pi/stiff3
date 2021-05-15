# stiff3

`stiff3` is a Fortran subprogram for solving stiff autonomous systems of ordinary differential equations (ODE's) using a semi-implicit Runge-Kutta method with three steps (SIRK3). The `stiff3` code was originally published in the work:

> Villadsen, J., & Michelsen, M. L. (1978). *Solution of differential equation models by polynomial approximation*. Prentice-Hall, Inc.

This repository provides a minimally refactored version with a modern procedural interface.

## Usage


## Requirements

Minimal requirements include:
* a recent Fortran compiler
* BLAS and LAPACK libraries

## Method


## Contributing

Feel welcome to propose additional features and submit bug reports by opening a new issue.

Here are some contribution ideas:
- Support for banded or sparse Jacobian matrices
- Use BLAS kernels for vector operations
- Continuous (dense) output of variables
- Extend `stiff3` to non-autonomous systems of ODE's
- Advanced stepsize control settings
- Improve the repository and code documentation
- Write a tutorial on how to use `stiff3`