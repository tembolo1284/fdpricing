/**
 * solver.h - PDE solver interface
 * 
 * Implements various finite difference schemes:
 * - Explicit Euler: Fast but requires small timesteps (CFL condition)
 * - Implicit Euler: Unconditionally stable, requires linear solve
 * - Crank-Nicolson: Second-order accurate in time
 * - PSOR: Projected SOR for American options with free boundary
 */

#ifndef FDPRICING_INTERNAL_SOLVER_H
#define FDPRICING_INTERNAL_SOLVER_H

#include "fdpricing.h"
#include "internal/core/context.h"
#include "internal/core/grid.h"
#include "internal/models/model.h"
#include "internal/options/option.h"

/* Solver parameters structure (internal) */
struct fdp_solver_params_s {
    fdp_context_t* ctx;
    fdp_solver_method_t method;
    double theta;           /* Theta parameter: 0=explicit, 0.5=CN, 1.0=implicit */
    double tolerance;       /* Convergence tolerance for iterative methods */
    int max_iterations;     /* Max iterations for PSOR */
    double omega;           /* SOR relaxation parameter (1.0-2.0) */
};

/* Result structure (internal) */
struct fdp_result_s {
    fdp_context_t* ctx;
    
    /* Solution surface: V(S, t) stored as V[time_idx * n_space + space_idx] */
    double* surface;
    int n_space;
    int n_time;
    
    /* Grid for interpolation */
    const double* space_points;
    const double* time_points;
    
    /* Convergence info */
    int iterations;
    fdp_error_t error_code;
    
    /* Reference to grid (not owned) */
    const fdp_grid_t* grid;
};

/* ========================================================================
 * Solver Parameter Management
 * ======================================================================== */

/* Already have public API, these are just helpers */

/* ========================================================================
 * 1D PDE Solvers
 * ======================================================================== */

/**
 * Solve 1D PDE using specified method
 * 
 * Solves backward in time from maturity to present
 * Returns result object with full solution surface
 */
fdp_result_t* fdp_solve_1d(
    fdp_context_t* ctx,
    const fdp_model_t* model,
    const fdp_option_t* option,
    const fdp_grid_t* grid,
    const fdp_solver_params_t* params
);

/**
 * Explicit Euler step (forward difference in time)
 * CFL condition: dt <= (dS)^2 / (sigma^2 * S^2)
 */
void fdp_solver_explicit_step(
    double* V_new,           /* Output: solution at t + dt */
    const double* V_old,     /* Input: solution at t */
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate
);

/**
 * Implicit Euler step (backward difference in time)
 * Unconditionally stable, requires tridiagonal solve
 */
void fdp_solver_implicit_step(
    double* V_new,
    const double* V_old,
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate
);

/**
 * Crank-Nicolson step (theta = 0.5)
 * Second-order accurate in time
 */
void fdp_solver_crank_nicolson_step(
    double* V_new,
    const double* V_old,
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate
);

/**
 * PSOR step for American options
 * Projected Successive Over-Relaxation
 */
void fdp_solver_psor_step(
    double* V_new,
    const double* V_old,
    const double* exercise_value,  /* Intrinsic value at each grid point */
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate,
    double omega,
    double tolerance,
    int max_iterations,
    int* iterations_taken
);

/* ========================================================================
 * Boundary Conditions
 * ======================================================================== */

/**
 * Apply boundary conditions to solution vector
 */
void fdp_apply_boundary_conditions(
    double* V,
    const fdp_grid_t* grid,
    const fdp_option_t* option,
    double t,
    double rate
);

/* ========================================================================
 * Terminal Condition (Payoff)
 * ======================================================================== */

/**
 * Set terminal condition V(S, T) = payoff(S)
 */
void fdp_set_terminal_condition(
    double* V,
    const fdp_grid_t* grid,
    const fdp_option_t* option
);

/* ========================================================================
 * CFL Condition Check
 * ======================================================================== */

/**
 * Check CFL stability condition for explicit method
 * Returns 1 if stable, 0 if unstable
 */
int fdp_check_cfl_condition(
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double dt
);

#endif /* FDPRICING_INTERNAL_SOLVER_H */
