/* tests/test_black_scholes_validation.c - Validate against Black-Scholes formula */

#include <criterion/criterion.h>
#include "fdpricing.h"
#include <math.h>
#include <stdio.h>

/* ========================================================================
 * Black-Scholes Analytical Formulas (for validation)
 * ======================================================================== */

static double norm_cdf(double x)
{
    /* Approximation of cumulative normal distribution */
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

static double black_scholes_call(
    double S, double K, double r, double q, double sigma, double T)
{
    double d1 = (log(S/K) + (r - q + 0.5*sigma*sigma)*T) / (sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    
    return S * exp(-q*T) * norm_cdf(d1) - K * exp(-r*T) * norm_cdf(d2);
}

static double black_scholes_put(
    double S, double K, double r, double q, double sigma, double T)
{
    double d1 = (log(S/K) + (r - q + 0.5*sigma*sigma)*T) / (sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    
    return K * exp(-r*T) * norm_cdf(-d2) - S * exp(-q*T) * norm_cdf(-d1);
}

/* ========================================================================
 * Validation Tests
 * ======================================================================== */

Test(black_scholes, atm_call_validation)
{
    double S = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20, T = 1.0;
    
    double analytical = black_scholes_call(S, K, r, q, sigma, T);
    double numerical = fdp_price_european_call(S, K, r, q, sigma, T, 200, 200);
    
    cr_assert_float_eq(numerical, analytical, 0.5,
                       "ATM call: FD=%.4f, BS=%.4f, diff=%.4f",
                       numerical, analytical, fabs(numerical - analytical));
}

Test(black_scholes, atm_put_validation)
{
    double S = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20, T = 1.0;
    
    double analytical = black_scholes_put(S, K, r, q, sigma, T);
    double numerical = fdp_price_european_put(S, K, r, q, sigma, T, 200, 200);
    
    cr_assert_float_eq(numerical, analytical, 0.5,
                       "ATM put: FD=%.4f, BS=%.4f, diff=%.4f",
                       numerical, analytical, fabs(numerical - analytical));
}

Test(black_scholes, put_call_parity)
{
    double S = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20, T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, q, sigma, T, 200, 200);
    double put = fdp_price_european_put(S, K, r, q, sigma, T, 200, 200);
    
    double parity_lhs = call - put;
    double parity_rhs = S - K * exp(-r * T);
    
    cr_assert_float_eq(parity_lhs, parity_rhs, 0.5,
                       "Put-call parity: C-P=%.4f, S-Ke^(-rT)=%.4f",
                       parity_lhs, parity_rhs);
}

Test(black_scholes, itm_call_vs_intrinsic)
{
    double S = 100.0, K = 80.0, r = 0.05, q = 0.0, sigma = 0.20, T = 1.0;
    
    double call = fdp_price_european_call(S, K, r, q, sigma, T, 200, 200);
    double intrinsic = S - K * exp(-r * T);
    
    cr_assert_gt(call, intrinsic,
                 "ITM call should be > discounted intrinsic: call=%.4f, intrinsic=%.4f",
                 call, intrinsic);
}

Test(black_scholes, detailed_validation)
{
    double S = 100.0, K = 100.0, r = 0.05, q = 0.0, sigma = 0.20, T = 1.0;
    
    double bs_call = black_scholes_call(S, K, r, q, sigma, T);
    double bs_put = black_scholes_put(S, K, r, q, sigma, T);
    
    double fd_call = fdp_price_european_call(S, K, r, q, sigma, T, 200, 200);
    double fd_put = fdp_price_european_put(S, K, r, q, sigma, T, 200, 200);
    
    printf("\n");
    printf("Black-Scholes vs Finite Difference:\n");
    printf("  BS Call:  %.6f\n", bs_call);
    printf("  FD Call:  %.6f\n", fd_call);
    printf("  Error:    %.6f (%.2f%%)\n", fabs(fd_call - bs_call), 
           100.0 * fabs(fd_call - bs_call) / bs_call);
    printf("\n");
    printf("  BS Put:   %.6f\n", bs_put);
    printf("  FD Put:   %.6f\n", fd_put);
    printf("  Error:    %.6f (%.2f%%)\n", fabs(fd_put - bs_put),
           100.0 * fabs(fd_put - bs_put) / bs_put);
    printf("\n");
    
    cr_assert_float_eq(fd_call, bs_call, 0.5, "Call prices should match within 0.5");
    cr_assert_float_eq(fd_put, bs_put, 0.5, "Put prices should match within 0.5");
}
