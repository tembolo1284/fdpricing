/**
 * solver1d.c - 1D PDE solver implementation
 * 
 * Implements finite difference schemes for 1D PDEs
 */

#include "fdpricing.h"
#include "internal/core/solver.h"
#include "internal/core/context.h"
#include "internal/core/grid.h"
#include "internal/models/model.h"
#include "internal/options/option.h"
#include "internal/numerics/tridiag.h"
#include "internal/utils/allocator.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* Diagnostic function to check for NaN/Inf in arrays */
static int check_array_valid(const double* arr, int n, const char* name) {
    for (int i = 0; i < n; i++) {
        if (isnan(arr[i]) || isinf(arr[i])) {
            printf("ERROR: %s[%d] = %f (invalid!)\n", name, i, arr[i]);
            return 0;
        }
        if (fabs(arr[i]) > 1e10) {
            printf("WARNING: %s[%d] = %f (very large!)\n", name, i, arr[i]);
        }
    }
    return 1;
}

/* ========================================================================
 * Terminal Condition Setup
 * ======================================================================== */

void fdp_set_terminal_condition(
    double* V,
    const fdp_grid_t* grid,
    const fdp_option_t* option)
{
    const double* S = grid->space_points;
    int n = grid->n_space;
    
    /* DEBUG: Print space grid endpoints */
    printf("DEBUG TERMINAL: n_space=%d\n", n);
    printf("DEBUG TERMINAL: S[0]=%.6f, S[%d]=%.6f\n", S[0], n-1, S[n-1]);
    printf("DEBUG TERMINAL: Strike K=%.6f, Type=%d (0=call, 1=put)\n", 
           option->strike, option->type);
    
    /* Set V(S, T) = payoff(S) */
    for (int i = 0; i < n; ++i) {
        V[i] = fdp_payoff_vanilla(option->type, S[i], option->strike);
    }
    
    /* DEBUG: Print payoffs at boundaries */
    printf("DEBUG TERMINAL: Payoff at S[0]=%.6f: V[0]=%.6f\n", S[0], V[0]);
    printf("DEBUG TERMINAL: Payoff at S[%d]=%.6f: V[%d]=%.6f\n", 
           n-1, S[n-1], n-1, V[n-1]);
    
    /* DEBUG: Print a few middle points */
    int mid = n / 2;
    printf("DEBUG TERMINAL: Payoff at S[%d]=%.6f: V[%d]=%.6f\n", 
           mid, S[mid], mid, V[mid]);
}

/* ========================================================================
 * Boundary Conditions
 * ======================================================================== */

void fdp_apply_boundary_conditions(
    double* V,
    const fdp_grid_t* grid,
    const fdp_option_t* option,
    double t,
    double rate)
{
    /* Lower boundary: S = S_min */
    V[0] = fdp_boundary_lower(option, t, rate);
    
    /* Upper boundary: S = S_max */
    V[grid->n_space - 1] = fdp_boundary_upper(option, 
                                               grid->space_points[grid->n_space - 1],
                                               t, rate);
}

/* ========================================================================
 * CFL Condition Check
 * ======================================================================== */

int fdp_check_cfl_condition(
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double dt)
{
    /* For stability of explicit method, we need:
     * dt <= (dS)^2 / (sigma^2 * S^2)
     * 
     * Check this at all grid points
     */
    
    const double* S = grid->space_points;
    int n = grid->n_space;
    
    for (int i = 1; i < n - 1; ++i) {
        double mu, sigma, r;
        model->vtable->get_coefficients_1d(model, S[i], 0.0, &mu, &sigma, &r);
        
        /* Get local grid spacing */
        double dS = (i < n - 1) ? (S[i + 1] - S[i]) : (S[i] - S[i - 1]);
        
        /* CFL condition: dt <= dS^2 / sigma^2 */
        if (sigma > 0.0) {
            double dt_max = dS * dS / (sigma * sigma);
            if (dt > dt_max) {
                return 0;  /* Unstable */
            }
        }
    }
    
    return 1;  /* Stable */
}

/* ========================================================================
 * Explicit Euler Step
 * ======================================================================== */

void fdp_solver_explicit_step(
    double* V_new,
    const double* V_old,
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate)
{
    (void)rate;
    const double* S = grid->space_points;
    int n = grid->n_space;
    
    /* Explicit scheme: V^{n+1}_i = V^n_i + dt * L(V^n)
     * where L is the spatial differential operator
     * 
     * Discretization:
     * dV/dt = (r - q)S * dV/dS + (1/2)*sigma^2*S^2 * d^2V/dS^2 - r*V
     * 
     * Using central differences:
     * dV/dS ≈ (V_{i+1} - V_{i-1}) / (S_{i+1} - S_{i-1})
     * d^2V/dS^2 ≈ (V_{i+1} - 2*V_i + V_{i-1}) / h^2
     */
    
    /* Boundaries are fixed by boundary conditions */
    V_new[0] = V_old[0];
    V_new[n - 1] = V_old[n - 1];
    
    /* Interior points */
    for (int i = 1; i < n - 1; ++i) {
        double mu, sigma, r;
        model->vtable->get_coefficients_1d(model, S[i], t, &mu, &sigma, &r);
        
        /* Grid spacing */
        double dS_plus = S[i + 1] - S[i];
        double dS_minus = S[i] - S[i - 1];
        // double dS_center = 0.5 * (dS_plus + dS_minus);
        
        /* First derivative coefficient */
        double dV_dS = (V_old[i + 1] - V_old[i - 1]) / (dS_plus + dS_minus);
        
        /* Second derivative coefficient */
        double d2V_dS2 = 2.0 * ((V_old[i + 1] - V_old[i]) / dS_plus -
                                (V_old[i] - V_old[i - 1]) / dS_minus) / 
                         (dS_plus + dS_minus);
        
        /* PDE: dV/dt = mu*dV/dS + (1/2)*sigma^2*d^2V/dS^2 - r*V */
        double dV_dt = mu * dV_dS + 0.5 * sigma * sigma * d2V_dS2 - r * V_old[i];
        
        /* Explicit Euler: V_new = V_old + dt * dV_dt */
        V_new[i] = V_old[i] + dt * dV_dt;
    }
}

/* ========================================================================
 * Implicit Euler Step
 * ======================================================================== */

void fdp_solver_implicit_step(
    double* V_new,
    const double* V_old,
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate)
{
    /* Implicit scheme: V^{n+1} - dt * L(V^{n+1}) = V^n
     * 
     * This leads to a tridiagonal system:
     * a_i * V_{i-1} + b_i * V_i + c_i * V_{i+1} = d_i
     */
    (void)rate;
    const double* S = grid->space_points;
    int n = grid->n_space;
    
    /* Allocate tridiagonal system */
    double* a = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);  /* Lower diagonal */
    double* b = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);  /* Main diagonal */
    double* c = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);  /* Upper diagonal */
    double* d = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);  /* RHS */
    
    if (!a || !b || !c || !d) {
        if (a) fdp_ctx_free(model->ctx, a);
        if (b) fdp_ctx_free(model->ctx, b);
        if (c) fdp_ctx_free(model->ctx, c);
        if (d) fdp_ctx_free(model->ctx, d);
        return;
    }
    
    /* Build tridiagonal system */
    for (int i = 0; i < n; ++i) {
        if (i == 0 || i == n - 1) {
            /* Boundary points: V_new[i] = V_old[i] (fixed by BC) */
            a[i] = 0.0;
            b[i] = 1.0;
            c[i] = 0.0;
            d[i] = V_old[i];
        } else {
            /* Interior points */
            double mu, sigma, r;
            model->vtable->get_coefficients_1d(model, S[i], t + dt, &mu, &sigma, &r);
            
            double dS_plus = S[i + 1] - S[i];
            double dS_minus = S[i] - S[i - 1];
            // double h_squared = dS_plus * dS_minus;
            
            /* Coefficients for implicit discretization */
            /* V_i - dt * [mu*(V_{i+1} - V_{i-1})/(dS+ + dS-) 
             *            + (1/2)*sigma^2*2*(V_{i+1}/dS+ - V_i*(1/dS+ + 1/dS-) + V_{i-1}/dS-)/(dS+ + dS-)
             *            - r*V_i] = V_old_i
             */
            
            double alpha = 0.5 * sigma * sigma;
            
            /* Lower diagonal coefficient (V_{i-1}) */
            a[i] = -dt * (-mu / (dS_plus + dS_minus) + 
                          alpha * 2.0 / (dS_minus * (dS_plus + dS_minus)));
            
            /* Upper diagonal coefficient (V_{i+1}) */
            c[i] = -dt * (mu / (dS_plus + dS_minus) + 
                          alpha * 2.0 / (dS_plus * (dS_plus + dS_minus)));
            
            /* Main diagonal coefficient (V_i) */
            b[i] = 1.0 - dt * (-alpha * 2.0 * (1.0 / dS_plus + 1.0 / dS_minus) / 
                               (dS_plus + dS_minus) - r);
            
            /* Right-hand side */
            d[i] = V_old[i];
        }
    }
    
    /* Solve tridiagonal system */
    fdp_solve_tridiagonal(model->ctx, a, b, c, d, V_new, n);
    
    /* Cleanup */
    fdp_ctx_free(model->ctx, a);
    fdp_ctx_free(model->ctx, b);
    fdp_ctx_free(model->ctx, c);
    fdp_ctx_free(model->ctx, d);
}

/* ========================================================================
 * Crank-Nicolson Step
 * ======================================================================== */

void fdp_solver_crank_nicolson_step(
    double* V_new,
    const double* V_old,
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate)
{
    /* Crank-Nicolson: Average of explicit and implicit
     * V^{n+1} - (dt/2)*L(V^{n+1}) = V^n + (dt/2)*L(V^n)
     * 
     * Theta-method with theta = 0.5
     */
    (void)rate;
    
    const double* S = grid->space_points;
    int n = grid->n_space;
    double theta = 0.5;
    
    /* DEBUG: Check inputs */
    printf("DEBUG CN: n=%d, dt=%.6f, theta=%.2f\n", n, dt, theta);
    printf("DEBUG CN: S[0]=%.2f, S[n-1]=%.2f\n", S[0], S[n-1]);
    printf("DEBUG CN: V_old[0]=%.6f, V_old[n-1]=%.6f\n", V_old[0], V_old[n-1]);
    
    /* Allocate tridiagonal system */
    double* a = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    double* b = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    double* c = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    double* d = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    
    if (!a || !b || !c || !d) {
        printf("DEBUG CN: ERROR - Memory allocation failed!\n");
        if (a) fdp_ctx_free(model->ctx, a);
        if (b) fdp_ctx_free(model->ctx, b);
        if (c) fdp_ctx_free(model->ctx, c);
        if (d) fdp_ctx_free(model->ctx, d);
        return;
    }
    
    /* Build system */
    for (int i = 0; i < n; ++i) {
        if (i == 0 || i == n - 1) {
            a[i] = 0.0;
            b[i] = 1.0;
            c[i] = 0.0;
            d[i] = V_old[i];
        } else {
            double mu, sigma, r;
            model->vtable->get_coefficients_1d(model, S[i], t + dt * theta, 
                                               &mu, &sigma, &r);
            
            /* DEBUG: Check model coefficients at mid point */
            if (i == n/2) {
                printf("DEBUG CN: At S[%d]=%.2f: mu=%.6f, sigma=%.6f, r=%.6f\n",
                       i, S[i], mu, sigma, r);
            }
            
            double dS_plus = S[i + 1] - S[i];
            double dS_minus = S[i] - S[i - 1];
            double alpha = 0.5 * sigma * sigma;
            
            /* DEBUG: Check grid spacing at mid point */
            if (i == n/2) {
                printf("DEBUG CN: dS_plus=%.6f, dS_minus=%.6f, alpha=%.6f\n",
                       dS_plus, dS_minus, alpha);
            }
            
            /* Implicit part coefficients (LHS) */
            double a_impl = -theta * dt * (-mu / (dS_plus + dS_minus) + 
                            alpha * 2.0 / (dS_minus * (dS_plus + dS_minus)));
            
            double c_impl = -theta * dt * (mu / (dS_plus + dS_minus) + 
                            alpha * 2.0 / (dS_plus * (dS_plus + dS_minus)));
            
            double b_impl = 1.0 - theta * dt * (-alpha * 2.0 * (1.0 / dS_plus + 1.0 / dS_minus) / 
                            (dS_plus + dS_minus) - r);
            
            /* Explicit part coefficients (RHS) */
            double a_expl = (1.0 - theta) * dt * (-mu / (dS_plus + dS_minus) + 
                            alpha * 2.0 / (dS_minus * (dS_plus + dS_minus)));
            
            double c_expl = (1.0 - theta) * dt * (mu / (dS_plus + dS_minus) + 
                            alpha * 2.0 / (dS_plus * (dS_plus + dS_minus)));
            
            double b_expl = 1.0 + (1.0 - theta) * dt * (-alpha * 2.0 * (1.0 / dS_plus + 1.0 / dS_minus) / 
                            (dS_plus + dS_minus) - r);
            
            /* DEBUG: Check matrix coefficients at mid point */
            if (i == n/2) {
                printf("DEBUG CN: LHS matrix row %d: a=%.6f, b=%.6f, c=%.6f\n", 
                       i, a_impl, b_impl, c_impl);
                printf("DEBUG CN: RHS coeffs: a_expl=%.6f, b_expl=%.6f, c_expl=%.6f\n",
                       a_expl, b_expl, c_expl);
            }
            
            /* Set system coefficients */
            a[i] = a_impl;
            b[i] = b_impl;
            c[i] = c_impl;
            
            /* RHS: Apply explicit part to V_old */
            d[i] = b_expl * V_old[i] + a_expl * V_old[i - 1] + c_expl * V_old[i + 1];
            
            /* DEBUG: Check RHS computation at mid point */
            if (i == n/2) {
                printf("DEBUG CN: d[%d] = %.6f*%.6f + %.6f*%.6f + %.6f*%.6f = %.6f\n",
                       i, b_expl, V_old[i], a_expl, V_old[i-1], c_expl, V_old[i+1], d[i]);
            }
        }
    }
    
    /* DEBUG: Check boundary conditions and a few RHS values */
    printf("DEBUG CN: d[0]=%.6f, d[%d]=%.6f, d[n-1]=%.6f\n", 
           d[0], n/2, d[n/2], d[n-1]);
    
    /* Solve */
    printf("DEBUG CN: About to solve tridiagonal system...\n");
    fdp_solve_tridiagonal(model->ctx, a, b, c, d, V_new, n);
    
    /* DEBUG: Check solution */
    printf("DEBUG CN: After solve: V_new[0]=%.6f, V_new[%d]=%.6f, V_new[n-1]=%.6f\n", 
           V_new[0], n/2, V_new[n/2], V_new[n-1]);
    
    /* Check for NaN/Inf */
    int has_invalid = 0;
    for (int i = 0; i < n; i++) {
        if (isnan(V_new[i]) || isinf(V_new[i])) {
            printf("DEBUG CN: ERROR - V_new[%d] = %f is invalid!\n", i, V_new[i]);
            has_invalid = 1;
            if (has_invalid > 5) break;  /* Don't spam too much */
        }
    }
    
    /* Cleanup */
    fdp_ctx_free(model->ctx, a);
    fdp_ctx_free(model->ctx, b);
    fdp_ctx_free(model->ctx, c);
    fdp_ctx_free(model->ctx, d);
}

/* ========================================================================
 * PSOR Step (Projected Successive Over-Relaxation)
 * ======================================================================== */

void fdp_solver_psor_step(
    double* V_new,
    const double* V_old,
    const double* exercise_value,
    const fdp_grid_t* grid,
    const fdp_model_t* model,
    double t,
    double dt,
    double rate,
    double omega,
    double tolerance,
    int max_iterations,
    int* iterations_taken)
{
    /* PSOR for American options: V >= exercise_value (free boundary)
     * 
     * Algorithm:
     * 1. Solve implicit step to get V_temp
     * 2. Project: V = max(V_temp, exercise_value)
     * 3. Iterate with SOR until convergence
     */
    
    const double* S = grid->space_points;
    int n = grid->n_space;
    
    /* Start with implicit solution */
    fdp_solver_implicit_step(V_new, V_old, grid, model, t, dt, rate);
    
    /* Project onto constraint */
    for (int i = 0; i < n; ++i) {
        if (V_new[i] < exercise_value[i]) {
            V_new[i] = exercise_value[i];
        }
    }
    
    /* Build system matrices for SOR iterations */
    double* a = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    double* b = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    double* c = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    double* rhs = FDP_CTX_ALLOC_ARRAY(model->ctx, double, n);
    
    if (!a || !b || !c || !rhs) {
        if (a) fdp_ctx_free(model->ctx, a);
        if (b) fdp_ctx_free(model->ctx, b);
        if (c) fdp_ctx_free(model->ctx, c);
        if (rhs) fdp_ctx_free(model->ctx, rhs);
        if (iterations_taken) *iterations_taken = 0;
        return;
    }
    
    /* Build coefficient matrix (same as implicit) */
    for (int i = 0; i < n; ++i) {
        if (i == 0 || i == n - 1) {
            a[i] = 0.0;
            b[i] = 1.0;
            c[i] = 0.0;
            rhs[i] = V_old[i];
        } else {
            double mu, sigma, r;
            model->vtable->get_coefficients_1d(model, S[i], t + dt, &mu, &sigma, &r);
            
            double dS_plus = S[i + 1] - S[i];
            double dS_minus = S[i] - S[i - 1];
            double alpha = 0.5 * sigma * sigma;
            
            a[i] = -dt * (-mu / (dS_plus + dS_minus) + 
                          alpha * 2.0 / (dS_minus * (dS_plus + dS_minus)));
            
            c[i] = -dt * (mu / (dS_plus + dS_minus) + 
                          alpha * 2.0 / (dS_plus * (dS_plus + dS_minus)));
            
            b[i] = 1.0 - dt * (-alpha * 2.0 * (1.0 / dS_plus + 1.0 / dS_minus) / 
                               (dS_plus + dS_minus) - r);
            
            rhs[i] = V_old[i];
        }
    }
    
    /* PSOR iterations */
    int iter;
    for (iter = 0; iter < max_iterations; ++iter) {
        double max_change = 0.0;
        
        for (int i = 1; i < n - 1; ++i) {
            /* Gauss-Seidel update */
            double V_gs = (rhs[i] - a[i] * V_new[i - 1] - c[i] * V_new[i + 1]) / b[i];
            
            /* SOR relaxation */
            double V_sor = (1.0 - omega) * V_new[i] + omega * V_gs;
            
            /* Projection onto constraint */
            double V_projected = (V_sor > exercise_value[i]) ? V_sor : exercise_value[i];
            
            /* Track convergence */
            double change = fabs(V_projected - V_new[i]);
            if (change > max_change) {
                max_change = change;
            }
            
            V_new[i] = V_projected;
        }
        
        /* Check convergence */
        if (max_change < tolerance) {
            break;
        }
    }
    
    if (iterations_taken) {
        *iterations_taken = iter + 1;
    }
    
    /* Cleanup */
    fdp_ctx_free(model->ctx, a);
    fdp_ctx_free(model->ctx, b);
    fdp_ctx_free(model->ctx, c);
    fdp_ctx_free(model->ctx, rhs);
}

/* ========================================================================
 * Main 1D Solver
 * ======================================================================== */

fdp_result_t* fdp_solve_1d(
    fdp_context_t* ctx,
    const fdp_model_t* model,
    const fdp_option_t* option,
    const fdp_grid_t* grid,
    const fdp_solver_params_t* params)
{
    if (!ctx || !model || !option || !grid || !params) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    int n_space = grid->n_space;
    int n_time = grid->n_time;
    double dt = grid->dt;
    
    /* Extract rate from model */
    double rate = 0.05;  /* Default fallback */
    if (model && model->vtable && model->vtable->get_coefficients_1d) {
        double mu_dummy, sigma_dummy;
        /* Get rate from model at spot price and t=0 */
        model->vtable->get_coefficients_1d(model, grid->spot, 0.0, 
                                           &mu_dummy, &sigma_dummy, &rate);
    }
    
    /* Check CFL condition for explicit method */
    if (params->method == FDP_SOLVER_EXPLICIT) {
        if (!fdp_check_cfl_condition(grid, model, dt)) {
            fdp_ctx_set_error(ctx, FDP_ERROR_STABILITY);
            return NULL;
        }
    }
    
    /* Allocate result */
    fdp_result_t* result = FDP_CTX_ALLOC_ARRAY(ctx, fdp_result_t, 1);
    if (!result) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Allocate solution surface */
    result->surface = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space * (n_time + 1));
    if (!result->surface) {
        fdp_ctx_free(ctx, result);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    result->ctx = ctx;
    result->n_space = n_space;
    result->n_time = n_time;
    result->space_points = grid->space_points;
    result->time_points = grid->time_points;
    result->grid = grid;
    result->iterations = 0;
    result->error_code = FDP_SUCCESS;
    
    /* Work arrays */
    double* V_old = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space);
    double* V_new = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space);
    double* exercise = NULL;
    
    if (option->style == FDP_STYLE_AMERICAN) {
        exercise = FDP_CTX_ALLOC_ARRAY(ctx, double, n_space);
    }
    
    if (!V_old || !V_new || (option->style == FDP_STYLE_AMERICAN && !exercise)) {
        if (V_old) fdp_ctx_free(ctx, V_old);
        if (V_new) fdp_ctx_free(ctx, V_new);
        if (exercise) fdp_ctx_free(ctx, exercise);
        fdp_result_free(result);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Set terminal condition */
    fdp_set_terminal_condition(V_old, grid, option);
    
    /* DEBUG: Check terminal condition */
    printf("DEBUG SOLVER: Checking terminal condition...\n");
    if (!check_array_valid(V_old, n_space, "V_terminal")) {
        printf("ERROR: Terminal condition has invalid values!\n");
    }
    printf("DEBUG SOLVER: Terminal condition range: [%.6f, %.6f]\n", 
           V_old[0], V_old[n_space-1]);
    
    /* Store terminal condition in surface */
    memcpy(result->surface + n_time * n_space, V_old, n_space * sizeof(double));
    
    /* Time-stepping backward from maturity to present */
    for (int n = n_time - 1; n >= 0; --n) {
        double t = grid->time_points[n];
        
        /* DEBUG: Print info for first time step */
        if (n == n_time - 1) {
            printf("DEBUG SOLVER: First time step, t=%.6f, dt=%.6f, rate=%.6f\n", t, dt, rate);
        }
        
        /* Compute exercise value for American options */
        if (option->style == FDP_STYLE_AMERICAN) {
            for (int i = 0; i < n_space; ++i) {
                exercise[i] = fdp_option_exercise_value(option, grid->space_points[i]);
            }
        }
        
        /* Apply PDE solver step */
        if (option->style == FDP_STYLE_AMERICAN && params->method == FDP_SOLVER_PSOR) {
            int iters;
            fdp_solver_psor_step(V_new, V_old, exercise, grid, model, t, dt,
                                rate,
                                params->omega, params->tolerance, 
                                params->max_iterations, &iters);
            result->iterations += iters;
        } else {
            switch (params->method) {
                case FDP_SOLVER_EXPLICIT:
                    fdp_solver_explicit_step(V_new, V_old, grid, model, t, dt, rate);
                    break;
                case FDP_SOLVER_IMPLICIT:
                    fdp_solver_implicit_step(V_new, V_old, grid, model, t, dt, rate);
                    break;
                case FDP_SOLVER_CRANK_NICOLSON:
                    fdp_solver_crank_nicolson_step(V_new, V_old, grid, model, t, dt, rate);
                    break;
                default:
                    fdp_solver_crank_nicolson_step(V_new, V_old, grid, model, t, dt, rate);
                    break;
            }
            
            /* Apply American constraint if needed */
            if (option->style == FDP_STYLE_AMERICAN) {
                for (int i = 0; i < n_space; ++i) {
                    if (V_new[i] < exercise[i]) {
                        V_new[i] = exercise[i];
                    }
                }
            }
        }
        
        /* DEBUG: Check after first step */
        if (n == n_time - 1) {
            printf("DEBUG SOLVER: After first step, checking V_new...\n");
            if (!check_array_valid(V_new, n_space, "V_new_step1")) {
                printf("ERROR: V_new has invalid values after first step!\n");
                /* Print first few and last few values */
                printf("V_new[0:4] = [%.6f, %.6f, %.6f, %.6f]\n", 
                       V_new[0], V_new[1], V_new[2], V_new[3]);
                printf("V_new[%d:%d] = [%.6f, %.6f, %.6f, %.6f]\n",
                       n_space-4, n_space-1,
                       V_new[n_space-4], V_new[n_space-3], 
                       V_new[n_space-2], V_new[n_space-1]);
            }
        }
        
        /* Apply boundary conditions */
        fdp_apply_boundary_conditions(V_new, grid, option, t, rate);
        
        /* Store in surface */
        memcpy(result->surface + n * n_space, V_new, n_space * sizeof(double));
        
        /* Swap arrays */
        double* temp = V_old;
        V_old = V_new;
        V_new = temp;
    }
    
    /* Cleanup */
    fdp_ctx_free(ctx, V_old);
    fdp_ctx_free(ctx, V_new);
    if (exercise) fdp_ctx_free(ctx, exercise);
    
    return result;
}
