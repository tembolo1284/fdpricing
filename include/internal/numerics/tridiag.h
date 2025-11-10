/**
 * tridiag.h - Tridiagonal matrix solver
 * 
 * Solves systems of the form: A*x = b
 * where A is tridiagonal:
 *   a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
 * 
 * Uses Thomas algorithm (specialized Gaussian elimination)
 * Complexity: O(n)
 */

#ifndef FDPRICING_INTERNAL_TRIDIAG_H
#define FDPRICING_INTERNAL_TRIDIAG_H

/**
 * Solve tridiagonal system: A*x = d
 * 
 * @param a Lower diagonal (size n-1, a[0] unused)
 * @param b Main diagonal (size n)
 * @param c Upper diagonal (size n-1, c[n-1] unused)
 * @param d Right-hand side (size n)
 * @param x Solution vector (size n, output)
 * @param n Size of system
 * 
 * Note: This modifies c and d arrays during computation
 */
void fdp_solve_tridiagonal(
    const double* a,
    const double* b,
    double* c,        /* Modified during solve */
    double* d,        /* Modified during solve */
    double* x,
    int n
);

/**
 * Solve tridiagonal system with separate work arrays
 * Does not modify input arrays
 */
void fdp_solve_tridiagonal_safe(
    const double* a,
    const double* b,
    const double* c,
    const double* d,
    double* x,
    double* work,    /* Work array, size 2*n */
    int n
);

#endif /* FDPRICING_INTERNAL_TRIDIAG_H */
