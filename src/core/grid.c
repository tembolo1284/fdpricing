/**
 * grid.c - Grid generation and manipulation
 */

#include "fdpricing.h"
#include "internal/core/grid.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"
#include <math.h>

/* ========================================================================
 * Grid Creation (Public API)
 * ======================================================================== */

fdp_grid_t* fdp_grid_new_1d(
    fdp_context_t* ctx,
    fdp_grid_type_t grid_type,
    int n_space,
    int n_time,
    double spot,
    double s_min,
    double s_max,
    double t_max)
{
    if (!ctx || n_space < 3 || n_time < 1) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    if (s_min < 0.0 || s_max <= s_min || spot <= 0.0 || t_max <= 0.0) {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    /* Allocate grid structure */
    fdp_grid_t* grid = FDP_CTX_ALLOC_ARRAY(ctx, fdp_grid_t, 1);
    if (!grid) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    grid->ctx = ctx;
    grid->dimension = FDP_GRID_DIM_1D;
    grid->grid_type = grid_type;
    grid->n_space = n_space;
    grid->n_time = n_time;
    grid->spot = spot;
    grid->s_min = s_min;
    grid->s_max = s_max;
    grid->t_max = t_max;
    grid->dt = t_max / (double)n_time;
    
    /* 2D fields not used */
    grid->vol_points = NULL;
    grid->n_vol = 0;
    grid->v_min = 0.0;
    grid->v_max = 0.0;
    grid->vol = 0.0;
    grid->dv = NULL;
    
    /* Allocate arrays */
    grid->space_points = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space);
    grid->time_points = FDP_CTX_ALLOC_ARRAY(ctx, double, n_time + 1);
    grid->ds = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space - 1);
    
    if (!grid->space_points || !grid->time_points || !grid->ds) {
        fdp_grid_free(grid);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Generate spatial grid based on type */
    switch (grid_type) {
        case FDP_GRID_UNIFORM:
            fdp_grid_create_uniform(grid->space_points, n_space, s_min, s_max);
            break;
            
        case FDP_GRID_LOG_UNIFORM:
            fdp_grid_create_log_uniform(grid->space_points, n_space, s_min, s_max);
            break;
            
        case FDP_GRID_SINH:
            /* Concentrate around spot price by default */
            grid->concentration_point = spot;
            grid->concentration_factor = 1.0;  /* Can be tuned */
            fdp_grid_create_sinh(grid->space_points, n_space, s_min, s_max,
                                spot, 1.0);
            break;
            
        case FDP_GRID_ADAPTIVE:
            /* For now, fall back to sinh - adaptive refinement TBD */
            grid->concentration_point = spot;
            grid->concentration_factor = 1.5;
            fdp_grid_create_sinh(grid->space_points, n_space, s_min, s_max,
                                spot, 1.5);
            break;
            
        default:
            fdp_grid_free(grid);
            fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
            return NULL;
    }
    
    /* Generate time grid (always uniform) */
    for (int i = 0; i <= n_time; ++i) {
        grid->time_points[i] = ((double)i / (double)n_time) * t_max;
    }
    
    /* Compute spacing arrays */
    fdp_grid_compute_spacing(grid->ds, grid->space_points, n_space);
    
    fdp_ctx_set_error(ctx, FDP_SUCCESS);
    return grid;
}

fdp_grid_t* fdp_grid_new_2d(
    fdp_context_t* ctx,
    int n_space,
    int n_vol,
    int n_time,
    double spot,
    double s_min,
    double s_max,
    double vol,
    double v_min,
    double v_max,
    double t_max)
{
    if (!ctx || n_space < 3 || n_vol < 3 || n_time < 1) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    if (s_min <= 0.0 || s_max <= s_min || v_min < 0.0 || v_max <= v_min) {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    /* Allocate grid structure */
    fdp_grid_t* grid = FDP_CTX_ALLOC_ARRAY(ctx, fdp_grid_t, 1);
    if (!grid) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    grid->ctx = ctx;
    grid->dimension = FDP_GRID_DIM_2D;
    grid->grid_type = FDP_GRID_UNIFORM;  /* 2D typically uses uniform */
    grid->n_space = n_space;
    grid->n_vol = n_vol;
    grid->n_time = n_time;
    grid->spot = spot;
    grid->s_min = s_min;
    grid->s_max = s_max;
    grid->vol = vol;
    grid->v_min = v_min;
    grid->v_max = v_max;
    grid->t_max = t_max;
    grid->dt = t_max / (double)n_time;
    
    /* Allocate arrays */
    grid->space_points = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space);
    grid->vol_points = FDP_CTX_ALLOC_ARRAY(ctx, double, n_vol);
    grid->time_points = FDP_CTX_ALLOC_ARRAY(ctx, double, n_time + 1);
    grid->ds = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space - 1);
    grid->dv = FDP_CTX_ALLOC_ARRAY(ctx, double, n_vol - 1);
    
    if (!grid->space_points || !grid->vol_points || !grid->time_points ||
        !grid->ds || !grid->dv) {
        fdp_grid_free(grid);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Generate grids - use log-uniform for space, uniform for vol */
    fdp_grid_create_log_uniform(grid->space_points, n_space, s_min, s_max);
    fdp_grid_create_uniform(grid->vol_points, n_vol, v_min, v_max);
    
    /* Generate time grid */
    for (int i = 0; i <= n_time; ++i) {
        grid->time_points[i] = ((double)i / (double)n_time) * t_max;
    }
    
    /* Compute spacing arrays */
    fdp_grid_compute_spacing(grid->ds, grid->space_points, n_space);
    fdp_grid_compute_spacing(grid->dv, grid->vol_points, n_vol);
    
    fdp_ctx_set_error(ctx, FDP_SUCCESS);
    return grid;
}

void fdp_grid_free(fdp_grid_t* grid)
{
    if (!grid) return;
    
    fdp_context_t* ctx = grid->ctx;
    
    if (grid->space_points) fdp_ctx_free(ctx, grid->space_points);
    if (grid->vol_points) fdp_ctx_free(ctx, grid->vol_points);
    if (grid->time_points) fdp_ctx_free(ctx, grid->time_points);
    if (grid->ds) fdp_ctx_free(ctx, grid->ds);
    if (grid->dv) fdp_ctx_free(ctx, grid->dv);
    
    fdp_ctx_free(ctx, grid);
}

/* ========================================================================
 * Grid Accessors
 * ======================================================================== */

int fdp_grid_get_n_space(const fdp_grid_t* grid)
{
    return grid ? grid->n_space : 0;
}

int fdp_grid_get_n_time(const fdp_grid_t* grid)
{
    return grid ? grid->n_time : 0;
}

double fdp_grid_get_spot(const fdp_grid_t* grid)
{
    return grid ? grid->spot : 0.0;
}

const double* fdp_grid_get_space_points(const fdp_grid_t* grid)
{
    return grid ? grid->space_points : NULL;
}

const double* fdp_grid_get_time_points(const fdp_grid_t* grid)
{
    return grid ? grid->time_points : NULL;
}

/* ========================================================================
 * Grid Generation Functions
 * ======================================================================== */

void fdp_grid_create_uniform(
    double* points,
    int n,
    double min,
    double max)
{
    if (!points || n < 2) return;
    
    double dx = (max - min) / (double)(n - 1);
    
    for (int i = 0; i < n; ++i) {
        points[i] = min + (double)i * dx;
    }
}

void fdp_grid_create_log_uniform(
    double* points,
    int n,
    double min,
    double max)
{
    if (!points || n < 2 || min <= 0.0) return;
    
    double log_min = log(min);
    double log_max = log(max);
    double d_log = (log_max - log_min) / (double)(n - 1);
    
    for (int i = 0; i < n; ++i) {
        points[i] = exp(log_min + (double)i * d_log);
    }
}

void fdp_grid_create_sinh(
    double* points,
    int n,
    double min,
    double max,
    double concentration_point,
    double concentration_factor)
{
    if (!points || n < 2) return;
    
    /* Map [min, max] to [-1, 1] with concentration_point at 0 */
    double c = concentration_point;
    double d = (max - min) / 2.0;
    double a = concentration_factor;  /* Controls concentration strength */
    
    double sinh_a = sinh(a);
    
    for (int i = 0; i < n; ++i) {
        /* Map index to [-1, 1] */
        double xi = -1.0 + 2.0 * (double)i / (double)(n - 1);
        
        /* Apply sinh transform */
        double eta = sinh(a * xi) / sinh_a;
        
        /* Map to [min, max] with center at concentration_point */
        points[i] = c + d * eta;
        
        /* Clamp to bounds (in case of numerical issues) */
        if (points[i] < min) points[i] = min;
        if (points[i] > max) points[i] = max;
    }
}

void fdp_grid_compute_spacing(
    double* ds,
    const double* points,
    int n)
{
    if (!ds || !points || n < 2) return;
    
    for (int i = 0; i < n - 1; ++i) {
        ds[i] = points[i + 1] - points[i];
    }
}

/* ========================================================================
 * Grid Utilities
 * ======================================================================== */

int fdp_grid_find_index(
    const double* points,
    int n,
    double value)
{
    if (!points || n < 1) return -1;
    
    /* Handle edge cases */
    if (value <= points[0]) return 0;
    if (value >= points[n - 1]) return n - 1;
    
    /* Binary search */
    int left = 0;
    int right = n - 1;
    
    while (right - left > 1) {
        int mid = (left + right) / 2;
        if (points[mid] <= value) {
            left = mid;
        } else {
            right = mid;
        }
    }
    
    /* Return closest index */
    if (fabs(points[left] - value) < fabs(points[right] - value)) {
        return left;
    }
    return right;
}

double fdp_grid_interpolate(
    const double* points,
    const double* values,
    int n,
    double x)
{
    if (!points || !values || n < 1) return 0.0;
    
    /* Handle edge cases */
    if (x <= points[0]) return values[0];
    if (x >= points[n - 1]) return values[n - 1];
    
    /* Find bracketing indices */
    int i = 0;
    while (i < n - 1 && points[i + 1] < x) {
        ++i;
    }
    
    /* Linear interpolation */
    double x0 = points[i];
    double x1 = points[i + 1];
    double y0 = values[i];
    double y1 = values[i + 1];
    
    double t = (x - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}
