/* src/numerics/tridiag.c - Tridiagonal solver (Thomas algorithm) */

#include "fdpricing.h"
#include "internal/numerics/tridiag.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"
#include <stdlib.h>

/* ========================================================================
 * Thomas Algorithm for Tridiagonal Systems
 * ======================================================================== */

void fdp_solve_tridiagonal(
    fdp_context_t* ctx,
    const double* a,    /* Lower diagonal [1..n-1] */
    const double* b,    /* Main diagonal [0..n-1] */
    const double* c,    /* Upper diagonal [0..n-2] */
    const double* d,    /* Right-hand side [0..n-1] */
    double* x,          /* Solution [0..n-1] */
    int n)
{
    if (n <= 0) return;
    
    /* Allocate working arrays - we need copies since we modify them */
    double* b_work = FDP_CTX_ALLOC_ARRAY(ctx, double, n);
    double* d_work = FDP_CTX_ALLOC_ARRAY(ctx, double, n);
    
    if (!b_work || !d_work) {
        if (b_work) FDP_CTX_FREE(ctx, b_work);
        if (d_work) FDP_CTX_FREE(ctx, d_work);
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return;
    }
    
    /* Copy b and d to working arrays */
    for (int i = 0; i < n; ++i) {
        b_work[i] = b[i];
        d_work[i] = d[i];
    }
    
    /* Forward elimination */
    for (int i = 1; i < n; ++i) {
        double m = a[i - 1] / b_work[i - 1];
        b_work[i] = b_work[i] - m * c[i - 1];
        d_work[i] = d_work[i] - m * d_work[i - 1];
    }
    
    /* Back substitution */
    x[n - 1] = d_work[n - 1] / b_work[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = (d_work[i] - c[i] * x[i + 1]) / b_work[i];
    }
    
    /* Clean up working arrays */
    FDP_CTX_FREE(ctx, b_work);
    FDP_CTX_FREE(ctx, d_work);
}

/* ========================================================================
 * Tridiagonal Solver with Pivoting (more stable)
 * ======================================================================== */

void fdp_solve_tridiagonal_pivoting(
    fdp_context_t* ctx,
    const double* a,
    const double* b,
    const double* c,
    const double* d,
    double* x,
    int n)
{
    /* For now, just call the standard version */
    /* TODO: Implement partial pivoting for better numerical stability */
    fdp_solve_tridiagonal(ctx, a, b, c, d, x, n);
}
