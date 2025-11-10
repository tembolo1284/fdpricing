/* src/options/barrier.c - Barrier option pricing (knock-in/knock-out) */

#include "fdpricing.h"
#include "internal/core/context.h"
#include "internal/core/grid.h"
#include "internal/core/solver.h"
#include "internal/models/model.h"
#include "internal/options/option.h"
#include <math.h>
#include <stdlib.h>

/* ========================================================================
 * Up-and-Out Call Option
 * ======================================================================== */

double fdp_price_barrier_option(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    fdp_option_type_t option_type,
    fdp_barrier_type_t barrier_type,
    double rebate,
    int n_space,
    int n_time)
{
    /* Create context */
    fdp_context_t* ctx = fdp_context_new();
    if (!ctx) return -1.0;
    
    /* Create grid - ensure barrier is within domain */
    double s_min = fmin(0.0, barrier * 0.5);
    double s_max = fmax(3.0 * spot, barrier * 1.5);
    
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
    
    /* Create barrier option */
    fdp_option_t* option = fdp_option_new_barrier(
        ctx,
        option_type,
        strike,
        maturity,
        barrier_type,
        barrier,
        rebate
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

/* Convenience wrappers for specific barrier types */

double fdp_price_up_and_out_call(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    return fdp_price_barrier_option(
        spot, strike, barrier, rate, div_yield, vol, maturity,
        FDP_OPTION_CALL, FDP_BARRIER_UP_OUT, 0.0,
        n_space, n_time
    );
}

double fdp_price_down_and_out_put(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    return fdp_price_barrier_option(
        spot, strike, barrier, rate, div_yield, vol, maturity,
        FDP_OPTION_PUT, FDP_BARRIER_DOWN_OUT, 0.0,
        n_space, n_time
    );
}

double fdp_price_up_and_in_call(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    return fdp_price_barrier_option(
        spot, strike, barrier, rate, div_yield, vol, maturity,
        FDP_OPTION_CALL, FDP_BARRIER_UP_IN, 0.0,
        n_space, n_time
    );
}

double fdp_price_down_and_in_put(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time)
{
    return fdp_price_barrier_option(
        spot, strike, barrier, rate, div_yield, vol, maturity,
        FDP_OPTION_PUT, FDP_BARRIER_DOWN_IN, 0.0,
        n_space, n_time
    );
}
