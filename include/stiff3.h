#ifndef STIFF3_H_
#define STIFF3_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*stiff3_fun)(int neqn, const double y[neqn], double dy[neqn],
                           void *data);

typedef void (*stiff3_jac)(int neqn, const double y[neqn], double df[],
                           void *data);

typedef void (*stiff3_out)(double x, int neqn, const double y[neqn], int iha,
                           double qa, void *data);

/*
 * integrate a stiff system of autonomous ODEs
 */
void stiff3_c(int /* n */, stiff3_fun /* fun */, stiff3_jac /* dfun */,
              stiff3_out /* out */, int /* nprint */, double /* x0 */,
              double /* x1 */, double * /* h0 */, double /* eps */,
              double /* w */[], double /* y */[], void * /* udata */);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // STIFF3_H_