/* src/options/american.c - American option pricing with early exercise */

#include "fdpricing.h"
#include "internal/core/context.h"
#include "internal/core/grid.h"
#include "internal/core/solver.h"
#include "internal/models/model.h"
#include "internal/options/option.h"
#include <math.h>
#include <stdlib.h>

/* ========================================================================
 * American Call Option
 * ======================================================================== */

double fdp_price_american_call(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    /* Create context */
    fdp_context_t* ctx = fdp_context_new();
    if (!ctx) return -1.0;
    
    /* Create grid */
    double s_min = 0.0;
    double s_max = 3.0 * fmax(spot, strike);
    
    fdp_grid_t* grid = fdp_grid_new_1d(
        ctx,
        FDP_GRID_UNIFORM,
        n_space,
        n_time,
        spot,
        s_min,
        s_max,
        maturity
    );
    
    if (!grid) {
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Create GBM model */
    fdp_model_t* model = fdp_model_new_gbm(ctx, rate, div_yield, vol);
    if (!model) {
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Create American call option */
    fdp_option_t* option = fdp_option_new_vanilla(
        ctx,
        FDP_OPTION_CALL,
        FDP_STYLE_AMERICAN,
        strike,
        maturity
    );
    
    if (!option) {
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Create solver parameters - use PSOR for American options */
    fdp_solver_params_t* params = fdp_solver_params_new(ctx);
    if (!params) {
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    fdp_solver_params_set_method(params, FDP_SOLVER_PSOR);
    fdp_solver_params_set_omega(params, 1.2);        /* SOR parameter */
    fdp_solver_params_set_tolerance(params, 1e-6);
    fdp_solver_params_set_max_iterations(params, 100);
    
    /* Solve PDE */
    fdp_result_t* result = fdp_solve_pde(ctx, model, option, grid, params);
    
    if (!result) {
        fdp_solver_params_free(params);
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Extract price */
    double price = fdp_result_get_price(result, spot);
    
    /* Cleanup */
    fdp_result_free(result);
    fdp_solver_params_free(params);
    fdp_option_free(option);
    fdp_model_free(model);
    fdp_grid_free(grid);
    fdp_context_free(ctx);
    
    return price;
}

/* ========================================================================
 * American Put Option
 * ======================================================================== */

double fdp_price_american_put(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    /* Create context */
    fdp_context_t* ctx = fdp_context_new();
    if (!ctx) return -1.0;
    
    /* Create grid */
    double s_min = 0.0;
    double s_max = 3.0 * fmax(spot, strike);
    
    fdp_grid_t* grid = fdp_grid_new_1d(
        ctx,
        FDP_GRID_UNIFORM,
        n_space,
        n_time,
        spot,
        s_min,
        s_max,
        maturity
    );
    
    if (!grid) {
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Create GBM model */
    fdp_model_t* model = fdp_model_new_gbm(ctx, rate, div_yield, vol);
    if (!model) {
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Create American put option */
    fdp_option_t* option = fdp_option_new_vanilla(
        ctx,
        FDP_OPTION_PUT,
        FDP_STYLE_AMERICAN,
        strike,
        maturity
    );
    
    if (!option) {
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Create solver parameters */
    fdp_solver_params_t* params = fdp_solver_params_new(ctx);
    if (!params) {
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    fdp_solver_params_set_method(params, FDP_SOLVER_PSOR);
    fdp_solver_params_set_omega(params, 1.2);
    fdp_solver_params_set_tolerance(params, 1e-6);
    fdp_solver_params_set_max_iterations(params, 100);
    
    /* Solve PDE */
    fdp_result_t* result = fdp_solve_pde(ctx, model, option, grid, params);
    
    if (!result) {
        fdp_solver_params_free(params);
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    
    /* Extract price */
    double price = fdp_result_get_price(result, spot);
    
    /* Cleanup */
    fdp_result_free(result);
    fdp_solver_params_free(params);
    fdp_option_free(option);
    fdp_model_free(model);
    fdp_grid_free(grid);
    fdp_context_free(ctx);
    
    return price;
}
