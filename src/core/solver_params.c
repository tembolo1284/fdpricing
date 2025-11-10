/**
 * solver_params.c - Solver parameter management
 */

#include "fdpricing.h"
#include "internal/core/solver.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"

fdp_solver_params_t* fdp_solver_params_new(fdp_context_t* ctx)
{
    if (!ctx) return NULL;
    
    fdp_solver_params_t* params = FDP_CTX_ALLOC_ARRAY(ctx, fdp_solver_params_t, 1);
    if (!params) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Set default values */
    params->ctx = ctx;
    params->method = FDP_SOLVER_CRANK_NICOLSON;  /* Good default */
    params->theta = 0.5;                          /* Crank-Nicolson */
    params->tolerance = 1e-6;                     /* Convergence tolerance */
    params->max_iterations = 100;                 /* For PSOR */
    params->omega = 1.2;                          /* SOR relaxation */
    
    return params;
}

void fdp_solver_params_free(fdp_solver_params_t* params)
{
    if (!params) return;
    fdp_ctx_free(params->ctx, params);
}

void fdp_solver_params_set_method(
    fdp_solver_params_t* params,
    fdp_solver_method_t method)
{
    if (!params) return;
    
    params->method = method;
    
    /* Set theta based on method */
    switch (method) {
        case FDP_SOLVER_EXPLICIT:
            params->theta = 0.0;
            break;
        case FDP_SOLVER_IMPLICIT:
            params->theta = 1.0;
            break;
        case FDP_SOLVER_CRANK_NICOLSON:
            params->theta = 0.5;
            break;
        case FDP_SOLVER_PSOR:
            params->theta = 1.0;  /* PSOR is implicit-based */
            break;
    }
}

void fdp_solver_params_set_theta(
    fdp_solver_params_t* params,
    double theta)
{
    if (!params) return;
    
    /* Clamp to [0, 1] */
    if (theta < 0.0) theta = 0.0;
    if (theta > 1.0) theta = 1.0;
    
    params->theta = theta;
}

void fdp_solver_params_set_tolerance(
    fdp_solver_params_t* params,
    double tolerance)
{
    if (params && tolerance > 0.0) {
        params->tolerance = tolerance;
    }
}

void fdp_solver_params_set_max_iterations(
    fdp_solver_params_t* params,
    int max_iter)
{
    if (params && max_iter > 0) {
        params->max_iterations = max_iter;
    }
}

void fdp_solver_params_set_omega(
    fdp_solver_params_t* params,
    double omega)
{
    if (!params) return;
    
    /* SOR omega should be in (1, 2) for over-relaxation */
    /* Clamp to reasonable range */
    if (omega < 0.5) omega = 0.5;
    if (omega > 2.0) omega = 2.0;
    
    params->omega = omega;
}
