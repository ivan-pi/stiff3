"""STIFF3 solver of Villadsen & Michelsen"""

import numpy as np

from collections import namedtuple
from _stiff3 import stiff3_step

def oderks3_native(fun,jac,t0,t1,y0,eps,w,h0):
    """Variant with adaptive stepping"""

    t_arr = [t0]
    yarr = [ytmp]
    iha_arr = [0.0]
    qa_arr = [1.0]

    # Step counter and output callback
    k = 0
    def out(t,y,iha,qa):
        t_arr.append(t)
        yarr.append(y)
        iha_arr.append(iha)
        qa_arr.append(qa)
        k += 1

    # Integrate from x0 to x1, with output at every step
    output_at_end = False
    x0, x1 = t0, t1

    h = stiff3_step(
        fun, jac,
        out,
        output_at_end,
        t0, 1
        h,
        eps, w,
        ytmp
    )

    yarr = np.asarray(yarr)
    iha_arr = np.asarray(iha_arr)
    qa_arr = np.asarray(qa_arr)

    # Integration complete
    return tspan, yarr, iha_arr, qa_arr


def oderks3_fixed(fun,jac,tspan,y0,eps,w,h0):
    """Variant with output only at the knots"""

    y0 = np.asarray(y0, dtype=np.float64, order='C')
    ytmp = y0.copy()

    yarr = np.empty((len(y0),len(tspan)), dtype=np.float64, order='F')
    iha_arr = np.empty(len(tspan),dtype=np.int32)
    qa_arr = np.empty(len(tspan),dtype=np.float64)

    # Step counter and output callback
    k = 0
    def out(t,y,iha,qa):
        yarr[:,k] = y
        iha_arr[k] = iha
        qa_arr[k] = qa
        k += 1

    output_at_end = True

    for i in range(1,len(tspan)):
        x0, x1 = tspan[i-1:i]
        h = stiff3_step(
            fun, jac,
            out,
            output_at_end,
            x0, x1
            h,
            eps, w,
            ytmp
        )

    # Integration complete
    return tspan, yarr, iha_arr, qa_arr


def odestiff3(fun,jac,tspan,y0,*,eps,w,init_step=None):
    """Adaptive solver for stiff system of autonmous ODEs

    The method is based on a semi-implicit Runge-Kutta method
    of third order (SIRK3). The step-size selection is based on an
    adaptive bisection algorithm. The solver requires an accurate
    Jacobian matrix; numerical finite differences aren't suitable for
    this algorithm.

    Parameter
    ---------
    fun : callable
        Right-hand side function fun(t,y)
    jac : callable
        Jacobian matrix jac(t,y)
    tspan : array_like
        The time-knots where we'd like to evaluate the system
        of ODEs.
    y0 : array_like
        Initial condition vector.
    eps : float
        relative tolerance
    w : array_like
        absolute tolerance
    init_step : float, optional
        The size of the first step. If absent take a fraction of
        the interval length.
    """

    # FIXME: smarter way to choose initial step-size
    if init_step:
        h = float(init_step)
    else:
        h = 0.001*(tspan[-1] - tspan[0])

    def _fun(y,dy):
        dy[:] = fun(y):

    def _jac(y,df):
        df[:,:] = jac(y)

    y0 = np.asarray(y0, dtype=np.float64, order='C')
    ytmp = y0.copy()

    tspan = np.asarray(tspan, dtype=np.float64)
    if len(tspan) < 2:
        raise ValueError("tspan should have at least two points")
    # FIXME: check tspan is monotonic

    if eps < 0:
        raise ValueError("eps must be positive.")


    w = np.asarray(w, dtype=np.float64, order='C')
    if w.shape != y0.shape:
        raise ValueError("y0 and w must have same shape.")
    if np.any(w <= 0.0):
        raise ValueError("elements of w must be positiove.")

    if len(tspan) == 2:
        t0, t1 = tspan
        tout, yout, rest = oderks3_native(
            _fun, _jac,
            t0, t1,
            y0,
            eps, w,
            h
        )
    else:
        tout, yout, rest = oderks3_fixed(
            _fun, _jac,
            tspan,
            y0,
            eps, w,
            h
        )

    return tout, yout