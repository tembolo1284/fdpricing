/* src/options/asian.c - Asian option pricing (average price options) */

#include "fdpricing.h"
#include "internal/core/context.h"
#include "internal/core/grid.h"
#include "internal/core/solver.h"
#include "internal/models/model.h"
#include "internal/options/option.h"
#include <math.h>
#include <stdlib.h>

/* ========================================================================
 * Asian Call Option (Arithmetic Average)
 * ======================================================================== */

double fdp_price_asian_call(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_averaging_points,
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
    
    /* Create Asian call option */
    fdp_option_t* option = fdp_option_new_asian(
        ctx,
        FDP_OPTION_CALL,
        strike,
        maturity,
        n_averaging_points
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
    
    fdp_solver_params_set_method(params, FDP_SOLVER_CRANK_NICOLSON);
    
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
 * Asian Put Option (Arithmetic Average)
 * ======================================================================== */

double fdp_price_asian_put(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_averaging_points,
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
    
    /* Create Asian put option */
    fdp_option_t* option = fdp_option_new_asian(
        ctx,
        FDP_OPTION_PUT,
        strike,
        maturity,
        n_averaging_points
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
    
    fdp_solver_params_set_method(params, FDP_SOLVER_CRANK_NICOLSON);
    
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
