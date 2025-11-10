/**
 * tridiag.c - Tridiagonal matrix solver (Thomas algorithm)
 */

#include "internal/numerics/tridiag.h"

void fdp_solve_tridiagonal(
    const double* a,
    const double* b,
    double* c,
    double* d,
    double* x,
    int n)
{
    if (n <= 0) return;
    
    /* Forward sweep: eliminate lower diagonal */
    for (int i = 1; i < n; ++i) {
        double m = a[i] / b[i - 1];
        b[i] = b[i] - m * c[i - 1];
        d[i] = d[i] - m * d[i - 1];
    }
    
    /* Back substitution */
    x[n - 1] = d[n - 1] / b[n - 1];
    
    for (int i = n - 2; i >= 0; --i) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
}

void fdp_solve_tridiagonal_safe(
    const double* a,
    const double* b,
    const double* c,
    const double* d,
    double* x,
    double* work,
    int n)
{
    if (n <= 0) return;
    
    /* Use work array for modified diagonals */
    double* b_work = work;         /* First n elements */
    double* d_work = work + n;     /* Next n elements */
    
    /* Copy b and d to work arrays */
    for (int i = 0; i < n; ++i) {
        b_work[i] = b[i];
        d_work[i] = d[i];
    }
    
    /* Forward sweep */
    for (int i = 1; i < n; ++i) {
        double m = a[i] / b_work[i - 1];
        b_work[i] = b_work[i] - m * c[i - 1];
        d_work[i] = d_work[i] - m * d_work[i - 1];
    }
    
    /* Back substitution */
    x[n - 1] = d_work[n - 1] / b_work[n - 1];
    
    for (int i = n - 2; i >= 0; --i) {
        x[i] = (d_work[i] - c[i] * x[i + 1]) / b_work[i];
    }
}
