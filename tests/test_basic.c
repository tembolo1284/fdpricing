/* tests/test_basic.c - Basic functionality tests using Criterion */

#include <criterion/criterion.h>
#include <criterion/redirect.h>
#include "fdpricing.h"
#include <math.h>

/* ========================================================================
 * Test Suite: Version and Compatibility
 * ======================================================================== */

Test(version, version_number_valid) {
    uint32_t version = fdp_get_version();
    cr_assert_gt(version, 0, "Version number should be greater than 0");
}

Test(version, library_compatible) {
    int compatible = fdp_is_compatible();
    cr_assert_neq(compatible, 0, "Library should be compatible");
}

/* ========================================================================
 * Test Suite: Context Management
 * ======================================================================== */

Test(context, create_and_destroy) {
    fdp_context_t* ctx = fdp_context_new();
    cr_assert_not_null(ctx, "Context creation should succeed");
    
    fdp_context_free(ctx);
    /* If we get here without crashing, cleanup succeeded */
    cr_assert(1);
}

Test(context, initial_error_state) {
    fdp_context_t* ctx = fdp_context_new();
    cr_assert_not_null(ctx);
    
    fdp_error_t error = fdp_get_last_error(ctx);
    cr_assert_eq(error, FDP_SUCCESS, "Initial error state should be SUCCESS");
    
    fdp_context_free(ctx);
}

Test(context, error_strings) {
    const char* success_str = fdp_get_error_string(FDP_SUCCESS);
    cr_assert_not_null(success_str, "Error string should not be NULL");
    cr_assert_str_eq(success_str, "Success", "Error string should be 'Success'");
    
    const char* alloc_str = fdp_get_error_string(FDP_ERROR_ALLOCATION);
    cr_assert_not_null(alloc_str);
    
    const char* invalid_str = fdp_get_error_string(FDP_ERROR_INVALID_PARAM);
    cr_assert_not_null(invalid_str);
}

/* ========================================================================
 * Test Suite: European Options - Basic Pricing
 * ======================================================================== */

Test(european_options, call_price_properties) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, q, sigma, T, 100, 100);
    
    cr_assert_gt(call, 0.0, "Call price should be positive");
    cr_assert_lt(call, S, "Call price should be less than spot price");
}

Test(european_options, put_price_properties) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double put = fdp_price_european_put(S, K, r, q, sigma, T, 100, 100);
    
    cr_assert_gt(put, 0.0, "Put price should be positive");
    cr_assert_lt(put, K, "Put price should be less than strike price");
}

Test(european_options, put_call_parity) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, q, sigma, T, 100, 100);
    double put = fdp_price_european_put(S, K, r, q, sigma, T, 100, 100);
    
    /* Put-call parity: C - P = S - K*exp(-r*T) */
    double parity_lhs = call - put;
    double parity_rhs = S - K * exp(-r * T);
    
    cr_assert_float_eq(parity_lhs, parity_rhs, 0.1,
                       "Put-call parity should hold: C - P = S - K*exp(-r*T)");
}

Test(european_options, atm_call_put_relationship) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, q, sigma, T, 100, 100);
    double put = fdp_price_european_put(S, K, r, q, sigma, T, 100, 100);
    
    /* For ATM options with no dividends, call should be slightly more expensive */
    cr_assert_gt(call, put, "ATM call should be worth more than ATM put");
}

/* ========================================================================
 * Test Suite: American Options
 * ======================================================================== */

Test(american_options, call_no_dividends) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double amer_call = fdp_price_american_call(S, K, r, q, sigma, T, 100, 100);
    double euro_call = fdp_price_european_call(S, K, r, q, sigma, T, 100, 100);
    
    cr_assert_gt(amer_call, 0.0, "American call price should be positive");
    cr_assert_float_eq(amer_call, euro_call, 0.1,
                       "American call should equal European call (no dividends)");
}

Test(american_options, put_early_exercise) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double amer_put = fdp_price_american_put(S, K, r, q, sigma, T, 100, 100);
    double euro_put = fdp_price_european_put(S, K, r, q, sigma, T, 100, 100);
    
    cr_assert_gt(amer_put, 0.0, "American put price should be positive");
    cr_assert_geq(amer_put, euro_put, 
                  "American put should be >= European put");
}

Test(american_options, put_itm_early_exercise_premium) {
    double S = 80.0;   /* ITM put */
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double amer_put = fdp_price_american_put(S, K, r, q, sigma, T, 150, 150);
    double euro_put = fdp_price_european_put(S, K, r, q, sigma, T, 150, 150);
    
    double premium = amer_put - euro_put;
    cr_assert_gt(premium, 0.1,
                 "ITM American put should have meaningful early exercise premium (>0.1)");
}

/* ========================================================================
 * Test Suite: Monotonicity Properties
 * ======================================================================== */

Test(monotonicity, call_increases_with_volatility) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double T = 1.0;
    
    double call_low_vol = fdp_price_european_call(S, K, r, q, 0.15, T, 100, 100);
    double call_high_vol = fdp_price_european_call(S, K, r, q, 0.25, T, 100, 100);
    
    cr_assert_gt(call_high_vol, call_low_vol,
                 "Call price should increase with volatility");
}

Test(monotonicity, call_increases_with_maturity) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    
    double call_short = fdp_price_european_call(S, K, r, q, sigma, 0.5, 100, 100);
    double call_long = fdp_price_european_call(S, K, r, q, sigma, 1.5, 100, 100);
    
    cr_assert_gt(call_long, call_short,
                 "Call price should increase with maturity");
}

Test(monotonicity, call_increases_with_spot) {
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double call_low_spot = fdp_price_european_call(90.0, K, r, q, sigma, T, 100, 100);
    double call_high_spot = fdp_price_european_call(110.0, K, r, q, sigma, T, 100, 100);
    
    cr_assert_gt(call_high_spot, call_low_spot,
                 "Call price should increase with spot price");
}

Test(monotonicity, put_increases_with_volatility) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double T = 1.0;
    
    double put_low_vol = fdp_price_european_put(S, K, r, q, 0.15, T, 100, 100);
    double put_high_vol = fdp_price_european_put(S, K, r, q, 0.25, T, 100, 100);
    
    cr_assert_gt(put_high_vol, put_low_vol,
                 "Put price should increase with volatility");
}

Test(monotonicity, put_decreases_with_spot) {
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double put_low_spot = fdp_price_european_put(90.0, K, r, q, sigma, T, 100, 100);
    double put_high_spot = fdp_price_european_put(110.0, K, r, q, sigma, T, 100, 100);
    
    cr_assert_lt(put_high_spot, put_low_spot,
                 "Put price should decrease with spot price");
}

/* ========================================================================
 * Test Suite: Boundary Cases
 * ======================================================================== */

Test(boundary, deep_itm_call_approaches_intrinsic) {
    double S = 150.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, q, sigma, T, 100, 100);
    double intrinsic = S - K * exp(-r * T);
    
    cr_assert_float_eq(call, intrinsic, 2.0,
                       "Deep ITM call should approach intrinsic value");
}

Test(boundary, deep_otm_call_near_zero) {
    double S = 50.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, q, sigma, T, 100, 100);
    
    cr_assert_lt(call, 1.0, "Deep OTM call should be near zero");
}

Test(boundary, deep_itm_put_approaches_intrinsic) {
    double S = 50.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double put = fdp_price_european_put(S, K, r, q, sigma, T, 100, 100);
    double intrinsic = K * exp(-r * T) - S;
    
    cr_assert_float_eq(put, intrinsic, 2.0,
                       "Deep ITM put should approach intrinsic value");
}

Test(boundary, deep_otm_put_near_zero) {
    double S = 150.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double put = fdp_price_european_put(S, K, r, q, sigma, T, 100, 100);
    
    cr_assert_lt(put, 1.0, "Deep OTM put should be near zero");
}

/* ========================================================================
 * Test Suite: Grid Convergence (Quick Check)
 * ======================================================================== */

Test(convergence, finer_grid_stability) {
    double S = 100.0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.0;
    double sigma = 0.20;
    double T = 1.0;
    
    double call_coarse = fdp_price_european_call(S, K, r, q, sigma, T, 50, 50);
    double call_fine = fdp_price_european_call(S, K, r, q, sigma, T, 200, 200);
    
    /* Prices should be similar (within 5%) */
    double rel_diff = fabs(call_fine - call_coarse) / call_coarse;
    cr_assert_lt(rel_diff, 0.05,
                 "Coarse and fine grid should produce similar results");
}
