dV/**
 * solver.c - Main solver dispatch
 */

#include "fdpricing.h"
#include "internal/core/solver.h"
#include "internal/core/context.h"
#include "internal/models/model.h"

fdp_result_t* fdp_solve_pde(
    fdp_context_t* ctx,
    const fdp_model_t* model,
    const fdp_option_t* option,
    const fdp_grid_t* grid,
    const fdp_solver_params_t* params)
{
    if (!ctx || !model || !option || !grid) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    /* Use default params if not provided */
    fdp_solver_params_t* default_params = NULL;
    if (!params) {
        default_params = fdp_solver_params_new(ctx);
        if (!default_params) {
            return NULL;
        }
        params = default_params;
    }
    
    /* Dispatch based on model dimensionality */
    fdp_result_t* result = NULL;
    
    if (fdp_model_is_1d(model)) {
        result = fdp_solve_1d(ctx, model, option, grid, params);
    } else if (fdp_model_is_2d(model)) {
        /* 2D solver not yet implemented */
        fdp_ctx_set_error(ctx, FDP_ERROR_NOT_IMPLEMENTED);
        result = NULL;
    } else {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        result = NULL;
    }
    
    /* Cleanup default params if allocated */
    if (default_params) {
        fdp_solver_params_free(default_params);
    }
    
    return result;
}
