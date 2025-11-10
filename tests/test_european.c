/* tests/test_european_criterion.c - European options tests with Criterion */

#include <criterion/criterion.h>
#include <criterion/parameterized.h>
#include "fdpricing.h"
#include <math.h>

/* ========================================================================
 * Basic European Option Tests
 * ======================================================================== */

Test(european_options, call_price_positive) {
    double call = fdp_price_european_call(100.0, 100.0, 0.05, 0.0, 0.20, 1.0, 100, 100);
    cr_assert_gt(call, 0.0, "Call price should be positive");
    cr_assert_lt(call, 100.0, "Call price should be less than spot");
}

Test(european_options, put_price_positive) {
    double put = fdp_price_european_put(100.0, 100.0, 0.05, 0.0, 0.20, 1.0, 100, 100);
    cr_assert_gt(put, 0.0, "Put price should be positive");
    cr_assert_lt(put, 100.0, "Put price should be less than strike");
}

Test(european_options, put_call_parity) {
    double S = 100.0, K = 100.0, r = 0.05, T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, 0.0, 0.20, T, 100, 100);
    double put = fdp_price_european_put(S, K, r, 0.0, 0.20, T, 100, 100);
    
    double parity_lhs = call - put;
    double parity_rhs = S - K * exp(-r * T);
    
    cr_assert_float_eq(parity_lhs, parity_rhs, 0.1,
                       "Put-call parity: C - P = S - K*exp(-r*T)");
}

/* ========================================================================
 * Parametric Tests - Multiple Strikes
 * ======================================================================== */

struct strike_test {
    double strike;
    double expected_min_call;
};

ParameterizedTestParameters(european_options, call_strikes) {
    static struct strike_test params[] = {
        {80.0, 20.0},   /* Deep ITM */
        {90.0, 10.0},   /* ITM */
        {100.0, 5.0},   /* ATM */
        {110.0, 2.0},   /* OTM */
        {120.0, 0.5},   /* Deep OTM */
    };
    
    size_t nb_params = sizeof(params) / sizeof(struct strike_test);
    return cr_make_param_array(struct strike_test, params, nb_params);
}

ParameterizedTest(struct strike_test *param, european_options, call_strikes) {
    double call = fdp_price_european_call(100.0, param->strike, 0.05, 0.0, 
                                          0.20, 1.0, 100, 100);
    
    cr_assert_gt(call, 0.0, "Call price for K=%.0f should be positive", 
                 param->strike);
    
    if (param->strike < 100.0) {
        /* ITM calls should have significant value */
        cr_assert_gt(call, param->expected_min_call,
                     "ITM call for K=%.0f should be > %.2f", 
                     param->strike, param->expected_min_call);
    }
}

/* ========================================================================
 * Monotonicity Tests
 * ======================================================================== */

Test(monotonicity, call_increases_with_volatility) {
    double call_low = fdp_price_european_call(100.0, 100.0, 0.05, 0.0, 
                                              0.15, 1.0, 100, 100);
    double call_high = fdp_price_european_call(100.0, 100.0, 0.05, 0.0, 
                                               0.25, 1.0, 100, 100);
    
    cr_assert_gt(call_high, call_low, 
                 "Call price should increase with volatility");
}

Test(monotonicity, call_increases_with_maturity) {
    double call_short = fdp_price_european_call(100.0, 100.0, 0.05, 0.0, 
                                                0.20, 0.5, 100, 100);
    double call_long = fdp_price_european_call(100.0, 100.0, 0.05, 0.0, 
                                               0.20, 1.5, 100, 100);
    
    cr_assert_gt(call_long, call_short,
                 "Call price should increase with maturity");
}

Test(monotonicity, call_increases_with_spot) {
    double call_low = fdp_price_european_call(90.0, 100.0, 0.05, 0.0, 
                                              0.20, 1.0, 100, 100);
    double call_high = fdp_price_european_call(110.0, 100.0, 0.05, 0.0, 
                                               0.20, 1.0, 100, 100);
    
    cr_assert_gt(call_high, call_low,
                 "Call price should increase with spot price");
}

/* ========================================================================
 * American Options Tests
 * ======================================================================== */

Test(american_options, put_early_exercise_premium) {
    double amer_put = fdp_price_american_put(80.0, 100.0, 0.05, 0.0, 
                                             0.20, 1.0, 150, 150);
    double euro_put = fdp_price_european_put(80.0, 100.0, 0.05, 0.0, 
                                             0.20, 1.0, 150, 150);
    
    cr_assert_gt(amer_put, euro_put,
                 "American put should have early exercise premium");
    
    double premium = amer_put - euro_put;
    cr_assert_gt(premium, 0.1,
                 "Early exercise premium should be meaningful for ITM puts");
}

Test(american_options, call_no_dividend_equals_european) {
    double amer_call = fdp_price_american_call(100.0, 100.0, 0.05, 0.0, 
                                               0.20, 1.0, 100, 100);
    double euro_call = fdp_price_european_call(100.0, 100.0, 0.05, 0.0, 
                                               0.20, 1.0, 100, 100);
    
    cr_assert_float_eq(amer_call, euro_call, 0.1,
                       "American call = European call (no dividends)");
}

/* ========================================================================
 * Boundary Conditions
 * ======================================================================== */

Test(boundary, deep_itm_call_approaches_intrinsic) {
    double S = 150.0, K = 100.0, r = 0.05, T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, 0.0, 0.20, T, 100, 100);
    double intrinsic = S - K * exp(-r * T);
    
    cr_assert_float_eq(call, intrinsic, 2.0,
                       "Deep ITM call should approach intrinsic value");
}

Test(boundary, deep_otm_call_near_zero) {
    double call = fdp_price_european_call(50.0, 100.0, 0.05, 0.0, 
                                          0.20, 1.0, 100, 100);
    
    cr_assert_lt(call, 1.0, "Deep OTM call should be near zero");
}
