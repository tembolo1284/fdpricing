/* tests/test_barrier.c - Barrier option tests using Criterion */

#include <criterion/criterion.h>
#include "fdpricing.h"
#include <math.h>

/* Small helper to build a GBM + grid + barrier option and price it via PDE */
static double price_barrier_option_up_out_call(
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
    fdp_context_t* ctx = fdp_context_new();
    if (!ctx) return -1.0;

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

    fdp_model_t* model = fdp_model_new_gbm(ctx, rate, div_yield, vol);
    if (!model) {
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }

    /* Up-and-out call, zero rebate */
    fdp_option_t* option = fdp_option_new_barrier(
        ctx,
        FDP_OPTION_CALL,
        strike,
        maturity,
        FDP_BARRIER_UP_OUT,
        barrier,
        0.0
    );
    if (!option) {
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }

    fdp_solver_params_t* params = fdp_solver_params_new(ctx);
    if (!params) {
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    fdp_solver_params_set_method(params, FDP_SOLVER_CRANK_NICOLSON);

    fdp_result_t* result = fdp_solve_pde(ctx, model, option, grid, params);
    if (!result) {
        fdp_solver_params_free(params);
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }

    double price = fdp_result_get_price(result, spot);

    fdp_result_free(result);
    fdp_solver_params_free(params);
    fdp_option_free(option);
    fdp_model_free(model);
    fdp_grid_free(grid);
    fdp_context_free(ctx);

    return price;
}

static double price_barrier_option_down_out_put(
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
    fdp_context_t* ctx = fdp_context_new();
    if (!ctx) return -1.0;

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

    fdp_model_t* model = fdp_model_new_gbm(ctx, rate, div_yield, vol);
    if (!model) {
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }

    /* Down-and-out put, zero rebate */
    fdp_option_t* option = fdp_option_new_barrier(
        ctx,
        FDP_OPTION_PUT,
        strike,
        maturity,
        FDP_BARRIER_DOWN_OUT,
        barrier,
        0.0
    );
    if (!option) {
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }

    fdp_solver_params_t* params = fdp_solver_params_new(ctx);
    if (!params) {
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }
    fdp_solver_params_set_method(params, FDP_SOLVER_CRANK_NICOLSON);

    fdp_result_t* result = fdp_solve_pde(ctx, model, option, grid, params);
    if (!result) {
        fdp_solver_params_free(params);
        fdp_option_free(option);
        fdp_model_free(model);
        fdp_grid_free(grid);
        fdp_context_free(ctx);
        return -1.0;
    }

    double price = fdp_result_get_price(result, spot);

    fdp_result_free(result);
    fdp_solver_params_free(params);
    fdp_option_free(option);
    fdp_model_free(model);
    fdp_grid_free(grid);
    fdp_context_free(ctx);

    return price;
}

/* ========================================================================
 * Test Suite: Barrier Options – Basic Relationships
 * ======================================================================== */

Test(barrier_options, up_and_out_call_less_than_vanilla) {
    double S     = 100.0;
    double K     = 100.0;
    double B     = 120.0;   /* Barrier above spot */
    double r     = 0.05;
    double q     = 0.0;
    double sigma = 0.20;
    double T     = 1.0;

    int n_space = 150;
    int n_time  = 150;

    double vanilla_call = fdp_price_european_call(S, K, r, q, sigma, T,
                                                  n_space, n_time);
    double barrier_call = fdp_price_up_and_out_call(S, K, B,
                                                           r, q, sigma, T,
                                                           n_space, n_time);

    cr_assert_gt(vanilla_call, 0.0, "Vanilla call should be positive");
    cr_assert_gt(barrier_call, 0.0, "Up-and-out call should be positive");
    cr_assert_leq(barrier_call, vanilla_call,
                  "Up-and-out call must be <= vanilla call (knock-out feature lowers value)");
}

Test(barrier_options, up_and_out_call_monotone_in_barrier) {
    double S     = 100.0;
    double K     = 100.0;
    double r     = 0.05;
    double q     = 0.0;
    double sigma = 0.20;
    double T     = 1.0;

    int n_space = 150;
    int n_time  = 150;

    /* Lower barrier => easier to knock out => lower value */
    double B_close = 110.0;
    double B_mid   = 130.0;
    double B_far   = 200.0;

    double price_close = fdp_price_up_and_out_call(
        S, K, r, q, sigma, T, B_close, 0.0, n_space, n_time);

    double price_mid = fdp_price__up_and_out_call(
        S, K, r, q, sigma, T, B_mid, 0.0, n_space, n_time);

    double price_far = fdp_price_up_and_out_call(
        S, K, r, q, sigma, T, B_far, 0.0, n_space, n_time);

    /* Allow tiny numerical noise – we only require *weak* monotonicity */
    double eps = 1e-3;

    cr_assert(price_mid >= price_close - eps,
              "Up-and-out call should not *decrease* as barrier moves further away "
              "(close=%.6f, mid=%.6f)",
              price_close, price_mid);

    cr_assert(price_far >= price_mid - eps,
              "Up-and-out call should not *decrease* as barrier moves further away "
              "(mid=%.6f, far=%.6f)",
              price_mid, price_far);
}

Test(barrier_options, down_and_out_put_less_than_vanilla) {
    double S     = 80.0;    /* Below strike */
    double K     = 100.0;
    double B     = 60.0;    /* Barrier below spot (down-and-out) */
    double r     = 0.05;
    double q     = 0.0;
    double sigma = 0.20;
    double T     = 1.0;

    int n_space = 150;
    int n_time  = 150;

    double vanilla_put = fdp_price_european_put(S, K, r, q, sigma, T,
                                                n_space, n_time);
    double barrier_put = fdp_price_down_out_put(S, K, B,
                                                           r, q, sigma, T,
                                                           n_space, n_time);

    cr_assert_gt(vanilla_put, 0.0, "Vanilla put should be positive");
    cr_assert_gt(barrier_put, 0.0, "Down-and-out put should be positive");
    cr_assert_leq(barrier_put, vanilla_put,
                  "Down-and-out put must be <= vanilla put (knock-out lowers value)");
}

