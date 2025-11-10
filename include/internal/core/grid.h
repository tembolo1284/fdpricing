/**
 * grid.h - Grid structure for finite difference discretization
 * 
 * Supports 1D and 2D grids with various spacing strategies:
 * - Uniform spacing
 * - Log-uniform (uniform in log-space)
 * - Sinh transformation (concentration around strike)
 * - Adaptive refinement
 */

#ifndef FDPRICING_INTERNAL_GRID_H
#define FDPRICING_INTERNAL_GRID_H

#include "fdpricing.h"
#include "internal/core/context.h"

/* Grid dimension */
typedef enum {
    FDP_GRID_DIM_1D = 1,
    FDP_GRID_DIM_2D = 2
} fdp_grid_dim_t;

/* Internal grid structure */
struct fdp_grid_s {
    fdp_context_t* ctx;           /* Parent context */
    fdp_grid_dim_t dimension;     /* 1D or 2D */
    fdp_grid_type_t grid_type;    /* Spacing type */
    
    /* 1D grid (always present) */
    double* space_points;         /* S grid (spot prices) */
    int n_space;                  /* Number of spatial points */
    double s_min;                 /* Minimum S */
    double s_max;                 /* Maximum S */
    double spot;                  /* Current spot price */
    
    /* 2D grid (only if dimension == 2D) */
    double* vol_points;           /* v grid (variance/volatility) */
    int n_vol;                    /* Number of vol points */
    double v_min;                 /* Minimum v */
    double v_max;                 /* Maximum v */
    double vol;                   /* Current vol/variance */
    
    /* Time grid */
    double* time_points;          /* t grid */
    int n_time;                   /* Number of time steps */
    double t_max;                 /* Maturity time */
    double dt;                    /* Time step size */
    
    /* Spatial steps (useful for FD schemes) */
    double* ds;                   /* ds[i] = s[i+1] - s[i] */
    double* dv;                   /* dv[i] = v[i+1] - v[i] (2D only) */
    
    /* Grid refinement info */
    double concentration_point;   /* Point where grid is concentrated (e.g., strike) */
    double concentration_factor;  /* How much to concentrate */
};

/* Grid creation helpers (internal) */

/**
 * Create uniform grid in linear space
 */
void fdp_grid_create_uniform(
    double* points,
    int n,
    double min,
    double max
);

/**
 * Create uniform grid in log-space
 * Useful for lognormal processes
 */
void fdp_grid_create_log_uniform(
    double* points,
    int n,
    double min,
    double max
);

/**
 * Create sinh-transformed grid
 * Concentrates points around a specific value (e.g., strike)
 * 
 * Transform: S = c + sinh(a*x)/sinh(a) where x ∈ [-1, 1]
 * This gives fine spacing near x=0 (which maps to concentration_point)
 */
void fdp_grid_create_sinh(
    double* points,
    int n,
    double min,
    double max,
    double concentration_point,
    double concentration_factor
);

/**
 * Compute spacing array: ds[i] = points[i+1] - points[i]
 * Array must be pre-allocated with size (n-1)
 */
void fdp_grid_compute_spacing(
    double* ds,
    const double* points,
    int n
);

/**
 * Find grid index closest to a given value using binary search
 */
int fdp_grid_find_index(
    const double* points,
    int n,
    double value
);

/**
 * Linear interpolation on grid
 */
double fdp_grid_interpolate(
    const double* points,
    const double* values,
    int n,
    double x
);

#endif /* FDPRICING_INTERNAL_GRID_H */
