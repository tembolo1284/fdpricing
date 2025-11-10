/**
 * result.c - Result object management
 */

#include "fdpricing.h"
#include "internal/core/solver.h"
#include "internal/core/context.h"
#include "internal/core/grid.h"
#include "internal/utils/allocator.h"
#include <math.h>

void fdp_result_free(fdp_result_t* result)
{
    if (!result) return;
    
    if (result->surface) {
        fdp_ctx_free(result->ctx, result->surface);
    }
    
    fdp_ctx_free(result->ctx, result);
}

double fdp_result_get_price(const fdp_result_t* result, double spot)
{
    if (!result || !result->surface) return 0.0;
    
    /* Price is at t=0 (last computed time step) */
    /* Interpolate on space grid */
    const double* V_now = result->surface;  /* First n_space elements */
    
    return fdp_grid_interpolate(result->space_points, V_now, 
                                 result->n_space, spot);
}

double fdp_result_get_delta(const fdp_result_t* result, double spot)
{
    if (!result || !result->surface) return 0.0;
    
    /* Delta = dV/dS, computed by finite difference */
    const double* V_now = result->surface;
    const double* S = result->space_points;
    int n = result->n_space;
    
    /* Find bracketing indices */
    int i = 0;
    while (i < n - 1 && S[i + 1] < spot) {
        ++i;
    }
    
    if (i >= n - 1) i = n - 2;
    
    /* Central difference approximation */
    double delta = (V_now[i + 1] - V_now[i]) / (S[i + 1] - S[i]);
    
    return delta;
}

double fdp_result_get_gamma(const fdp_result_t* result, double spot)
{
    if (!result || !result->surface) return 0.0;
    
    /* Gamma = d^2V/dS^2, computed by second-order finite difference */
    const double* V_now = result->surface;
    const double* S = result->space_points;
    int n = result->n_space;
    
    /* Find bracketing indices */
    int i = 0;
    while (i < n - 1 && S[i + 1] < spot) {
        ++i;
    }
    
    if (i == 0) i = 1;
    if (i >= n - 1) i = n - 2;
    
    /* Three-point formula for second derivative */
    double h1 = S[i] - S[i - 1];
    double h2 = S[i + 1] - S[i];
    
    /* Second derivative by finite difference */
    double gamma = 2.0 * ((V_now[i + 1] - V_now[i]) / h2 - 
                          (V_now[i] - V_now[i - 1]) / h1) / (h1 + h2);
    
    return gamma;
}

double fdp_result_get_theta(const fdp_result_t* result, double spot)
{
    if (!result || !result->surface) return 0.0;
    
    /* Theta = -dV/dt (negative of time derivative) */
    /* Approximate using first two time steps */
    
    if (result->n_time < 2) return 0.0;
    
    const double* V_0 = result->surface;  /* t = 0 */
    const double* V_1 = result->surface + result->n_space;  /* t = dt */
    
    double dt = result->time_points[1] - result->time_points[0];
    
    double V0 = fdp_grid_interpolate(result->space_points, V_0, 
                                      result->n_space, spot);
    double V1 = fdp_grid_interpolate(result->space_points, V_1, 
                                      result->n_space, spot);
    
    /* Theta = -(V1 - V0) / dt (negative because time decreases) */
    return -(V1 - V0) / dt;
}

int fdp_result_get_iterations(const fdp_result_t* result)
{
    return result ? result->iterations : 0;
}

fdp_error_t fdp_result_get_error_code(const fdp_result_t* result)
{
    return result ? result->error_code : FDP_ERROR_INVALID_PARAM;
}

const double* fdp_result_get_surface(const fdp_result_t* result)
{
    return result ? result->surface : NULL;
}
