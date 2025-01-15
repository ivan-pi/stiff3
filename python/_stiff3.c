// _stiff3.c

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "_stiff3.h"

#define PYERR_RET(errobj, message)                                             \
  {                                                                            \
    PyErr_SetString(errobj, message);                                          \
    return NULL;                                                               \
  }
#define PYERR(errobj, message)                                                 \
  { PyErr_SetString(errobj, message); }

// (!) Only for 1-d arrays (!)
#define CHECK_NPY_PROPS(ap, type, msg)                                         \
  do {                                                                         \
    if (!(PyArray_IS_C_CONTIGUOUS(ap)) || (PyArray_TYPE(ap) != type)) {        \
      PYERR(rkc_error, msg)                                                    \
    }                                                                          \
  } while (0)

struct rkc_callback_capsule {
  PyObject *fun;
  PyObject *jac;
  PyObject *out; // Use this for adaptive mode
};

static void fun_adaptor(int neqn, const double y[neqn], double dy[neqn],
                        void *data) {
  struct callbacks *c = data;
  assert(c->fun != NULL);

  //
  // 1. create NumPy ndarray instances
  //
  const npy_intp dims_y[1] = {neqn};
  PyArrayObject *ap_y = (PyArrayObject *)PyArray_SimpleNewFromData(
      1, dims_y, NPY_FLOAT64, (void *)y);
  if (!ap_y)
    PYERR(rkc_error, "Failed to create NumPy ndarray from pointer (*y)")
  PyArray_CLEARFLAGS(ap_y, NPY_ARRAY_WRITEABLE);

  const npy_intp dims_dy[1] = {neqn};
  PyArrayObject *ap_dy = (PyArrayObject *)PyArray_SimpleNewFromData(
      1, dims_dy, NPY_FLOAT64, (void *)dy);
  if (!ap_dy)
    PYERR(rkc_error, "Failed to create NumPy ndarray from pointer (*dy)")

  //
  // 2. invoke the callback function
  //    fun(y,dy) -> None
  //    (the function is expected to fill/mutate dy)
  //
  PyObject *res = PyObject_CallFunction(caps->fun, "OO", ap_y, ap_dy);
}

static void jac_adaptor(int neqn, const double y[neqn], double df[],
                        void *data) {
  const struct callbacks *c = data;

  // State vector
  const npy_intp dims[1] = {neqn};
  PyArrayObject *ap_y = (PyArrayObject *)PyArray_SimpleNewFromData(
      1, dims, NPY_FLOAT64, (void *)y);
  if (!ap_y)
    PYERR(rkc_error, "Failed to create ndarray");
  PyArray_CLEARFLAGS(ap_y, NPY_ARRAY_WRITEABLE);

  // Jacobian array
  const npy_intp dims_df[2] = {neqn, neqn};
  PyArrayObject *ap_df = (PyArrayObject *)PyArray_SimpleNewFromData(
      2, dims_df, NPY_FLOAT64, (void *)df);
  if (!ap_df)
    PYERR(rkc_error, "Failed to create ndarray from pointer (*jac)")

  // FIXME: tell Numpy it is Fortran-ordered array

  //
  // 2. invoke the callback
  //    jac(y,dfdy) -> float:
  //
  PyObject *res = PyObject_CallFunction(c->jac, "OO", ap_y, ap_df);
  if (res == NULL) {
    // FIXME, how to handle this?
    PYERR(rkc_error, "Failed to call spcrad callback.")
  };
}

static void out_adaptor(double x, int neqn, const double y[neqn], int iha, double qa,
                        void *data) {

  struct callbacks *c = data;

  const npy_intp dims[1] = {neqn};
  PyArrayObject *ap_y = (PyArrayObject *)PyArray_SimpleNewFromData(
      1, dims, NPY_FLOAT64, (void *)y);
  if (!ap_y)
    PYERR(rkc_error, "Failed to create ndarray");
  PyArray_CLEARFLAGS(ap_y, NPY_ARRAY_WRITEABLE);

  // Out function
  //      def fout(t,y,iha,qa) -> None
  PyObject *res = PyObject_CallFunction(c->out, "dOid", x, ap_y, iha, qa);
}

// neqn is inferred from y
static char stiff3_step_doc[] =
    "h = stiff3_step(fun,jac,out,nout,x0,x1,h0,eps,w,y)";

static PyObject *stiff3_step(PyObject *self, PyObject *args) {

  int imode; // (bool)
  double x0, x1, h, eps;

  // FIXME: add basic counters
  // int nfev, naccpt, nrejct;

  PyArrayObject *ap_y = NULL, *ap_w = NULL;
  struct callbacks c = {.fun = NULL, .jac = NULL, .out = NULL};

  // h = stiff3(fun,jac,out,imode,x0,x1,h0,eps,w,y)
  if (!(PyArg_ParseTuple(args, "OOOpddddO!O!:stiff3_step",
                         &c.fun,                            // O
                         &c.jac,                            // O
                         &c.out,                            // O
                         &imode,                            // p
                         &x0, &x1,                          // dd
                         &h,                                // d
                         &eps,                              // d
                         &PyArray_Type, (PyObject **)&ap_w, // O!
                         &PyArray_Type, (PyObject **)&ap_y, // O!
                         ))) {
    return NULL;
  }

  assert(ap_y != NULL);
  assert(ap_w != NULL);
  assert(c.fun != NULL);
  assert(c.jac != NULL);

  // Check contiguity and type correctness
  CHECK_NPY_PROPS(ap_y, NPY_FLOAT64,
                  "Argument (y) must be a contiguous array of type float64.");
  CHECK_NPY_PROPS(ap_w, NPY_FLOAT64,
                  "Argument (w) must be a contiguous array of type float64.");

  int neqn = PyArray_DIMS(ap_y)[0];

  // Unpack the raw values
  double *y = (double *)PyArray_DATA(ap_y);
  double *w = (double *)PyArray_DATA(ap_w);

  stiff3_c(neqn, &fun_adaptor, &jac_adaptor, &out_adaptor, (imode) ? INT_MAX : 1,
           x0, x1, &h, eps, w, y, (void *)c);

  PyObject *res = Py_BuildValue("d", h);
  if (res == NULL) {
    return NULL;
  }

  return res;
}

static struct PyMethodDef stiff3_methods[] = {
    {"stiff3_step", stiff3_step, METH_VARARGS, stiff3_step_doc},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_stiff3",                                           /* m_name */
    "semi-implicit RK solver for stiff autonomous ODEs", /* m_doc */
    -1,                                                  /* m_size */
    stiff3_methods,                                      /* m_methods */
    NULL,                                                /* m_slots */
    NULL,                                                /* m_traverse */
    NULL,                                                /* m_clear */
    NULL                                                 /* m_free */
};

PyMODINIT_FUNC PyInit__stiff3(void) {
  // Initialize NumPy
  import_array();
  PyObject *m = PyModule_Create(&module_def);
  if (m == NULL) {
    return NULL;
  }
  return m;
}
