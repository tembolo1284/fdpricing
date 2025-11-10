/* include/internal/numerics/tridiag.h - Tridiagonal solver */

#ifndef FDP_INTERNAL_TRIDIAG_H
#define FDP_INTERNAL_TRIDIAG_H

#include "fdpricing.h"

/* Forward declaration */
struct fdp_context_s;

/* ========================================================================
 * Tridiagonal System Solver (Thomas Algorithm)
 * ======================================================================== */

/**
 * Solve tridiagonal system: A*x = d
 * 
 * The matrix A has the form:
 *   b[0]  c[0]   0     0    ...   0
 *   a[0]  b[1]  c[1]   0    ...   0
 *    0    a[1]  b[2]  c[2]  ...   0
 *   ...   ...   ...   ...   ...  ...
 *    0     0     0    a[n-2] b[n-1]
 * 
 * Uses Thomas algorithm (specialized Gaussian elimination for tridiagonal systems)
 * Time complexity: O(n), Space complexity: O(n)
 * 
 * @param ctx   Context for memory allocation
 * @param a     Lower diagonal [a[0], ..., a[n-2]] (length n-1)
 * @param b     Main diagonal [b[0], ..., b[n-1]] (length n)
 * @param c     Upper diagonal [c[0], ..., c[n-2]] (length n-1)
 * @param d     Right-hand side [d[0], ..., d[n-1]] (length n)
 * @param x     Solution vector (output) [x[0], ..., x[n-1]] (length n)
 * @param n     Size of the system
 */
void fdp_solve_tridiagonal(
    fdp_context_t* ctx,
    const double* a,
    const double* b,
    const double* c,
    const double* d,
    double* x,
    int n);

/**
 * Solve tridiagonal system with partial pivoting (more stable)
 * 
 * Same interface as fdp_solve_tridiagonal, but uses pivoting for
 * better numerical stability when diagonal elements are small.
 */
void fdp_solve_tridiagonal_pivoting(
    fdp_context_t* ctx,
    const double* a,
    const double* b,
    const double* c,
    const double* d,
    double* x,
    int n);

#endif /* FDP_INTERNAL_TRIDIAG_H */
