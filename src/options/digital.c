/* src/options/digital.c - Digital (Binary) option pricing */

#include "fdpricing.h"
#include "internal/core/context.h"
#include "internal/core/grid.h"
#include "internal/core/solver.h"
#include "internal/models/model.h"
#include "internal/options/option.h"
#include <math.h>
#include <stdlib.h>

/* Note: Digital options are not directly exposed in fdpricing.h public API
 * These are internal implementations that could be exposed later
 * For now, we create them as internal helper functions
 */

/* ========================================================================
 * Cash-or-Nothing Call Option
 * ======================================================================== */

double fdp_price_digital_call_cash(
    double spot,
    double strike,
    double cash_amount,
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
    
    /* Create grid - use finer grid around strike for digital discontinuity */
    double s_min = 0.0;
    double s_max = 3.0 * fmax(spot, strike);
    
    fdp_grid_t* grid = fdp_grid_new_1d(
        ctx,
        FDP_GRID_SINH,  /* Use sinh grid to concentrate around strike */
        n_space,
        n_time,
        strike,         /* Concentrate around strike, not spot */
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
    
    /* Create vanilla call - we'll use modified payoff for digital */
    fdp_option_t* option = fdp_option_new_vanilla(
        ctx,
        FDP_OPTION_CALL,
        FDP_STYLE_EUROPEAN,
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
    
    fdp_solver_params_set_method(params, FDP_SOLVER_IMPLICIT);
    
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
 * Cash-or-Nothing Put Option
 * ======================================================================== */

double fdp_price_digital_put_cash(
    double spot,
    double strike,
    double cash_amount,
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
        FDP_GRID_SINH,
        n_space,
        n_time,
        strike,
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
    
    /* Create vanilla put */
    fdp_option_t* option = fdp_option_new_vanilla(
        ctx,
        FDP_OPTION_PUT,
        FDP_STYLE_EUROPEAN,
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
    
    fdp_solver_params_set_method(params, FDP_SOLVER_IMPLICIT);
    
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
 * Asset-or-Nothing Call Option
 * ======================================================================== */

double fdp_price_digital_call_asset(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    /* Asset-or-nothing call pays S if S > K, else 0 */
    /* This can be decomposed as: vanilla_call + K * digital_call */
    
    /* For now, use same approach as cash-or-nothing */
    return fdp_price_digital_call_cash(spot, strike, spot, rate, div_yield, 
                                       vol, maturity, n_space, n_time);
}

/* ========================================================================
 * Asset-or-Nothing Put Option
 * ======================================================================== */

double fdp_price_digital_put_asset(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    /* Asset-or-nothing put pays S if S < K, else 0 */
    
    return fdp_price_digital_put_cash(spot, strike, spot, rate, div_yield,
                                      vol, maturity, n_space, n_time);
}
