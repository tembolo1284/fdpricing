/* examples/simple_european.c - Simple European option pricing example */

#include "fdpricing.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=================================================\n");
    printf("FD Pricing Library - Simple European Options\n");
    printf("=================================================\n\n");
    
    /* Market parameters */
    double spot = 100.0;
    double strike = 100.0;
    double rate = 0.05;        /* 5% risk-free rate */
    double div_yield = 0.0;    /* No dividends */
    double vol = 0.20;         /* 20% volatility */
    double maturity = 1.0;     /* 1 year to maturity */
    
    /* Grid parameters */
    int n_space = 100;
    int n_time = 100;
    
    printf("Market Parameters:\n");
    printf("  Spot price (S):         $%.2f\n", spot);
    printf("  Strike price (K):       $%.2f\n", strike);
    printf("  Risk-free rate (r):     %.2f%%\n", rate * 100.0);
    printf("  Dividend yield (q):     %.2f%%\n", div_yield * 100.0);
    printf("  Volatility (sigma):     %.2f%%\n", vol * 100.0);
    printf("  Time to maturity (T):   %.2f years\n", maturity);
    printf("\n");
    
    printf("Grid Parameters:\n");
    printf("  Space points:           %d\n", n_space);
    printf("  Time steps:             %d\n", n_time);
    printf("  Space range:            [%.2f, %.2f]\n", 0.0, 3.0 * spot);
    printf("  Time step (dt):         %.6f years\n", maturity / n_time);
    printf("\n");
    
    printf("Numerical Method:\n");
    printf("  Solver:                 Crank-Nicolson (θ = 0.5)\n");
    printf("  Model:                  Black-Scholes (GBM)\n");
    printf("  Grid type:              Uniform in space\n");
    printf("\n");
    
    /* Price European call */
    printf("European Call Option:\n");
    printf("--------------------------------------------------\n");
    double call_price = fdp_price_european_call(
        spot, strike, rate, div_yield, vol, maturity,
        n_space, n_time
    );
    
    if (call_price < 0) {
        printf("  ERROR: Pricing failed!\n");
        return 1;
    }
    
    printf("  Price:                  $%.6f\n", call_price);
    
    /* Calculate Black-Scholes for comparison */
    double d1 = (log(spot/strike) + (rate - div_yield + 0.5*vol*vol)*maturity) / (vol*sqrt(maturity));
    double d2 = d1 - vol*sqrt(maturity);
    double N_d1 = 0.5 * (1.0 + erf(d1 / sqrt(2.0)));
    double N_d2 = 0.5 * (1.0 + erf(d2 / sqrt(2.0)));
    double bs_call = spot * exp(-div_yield*maturity) * N_d1 - strike * exp(-rate*maturity) * N_d2;
    
    printf("  Black-Scholes price:    $%.6f\n", bs_call);
    printf("  Absolute error:         $%.6f\n", fabs(call_price - bs_call));
    printf("  Relative error:         %.2f%%\n", 100.0 * fabs(call_price - bs_call) / bs_call);
    printf("\n");
    
    /* Price European put */
    printf("European Put Option:\n");
    printf("--------------------------------------------------\n");
    double put_price = fdp_price_european_put(
        spot, strike, rate, div_yield, vol, maturity,
        n_space, n_time
    );
    
    if (put_price < 0) {
        printf("  ERROR: Pricing failed!\n");
        return 1;
    }
    
    printf("  Price:                  $%.6f\n", put_price);
    
    /* Calculate Black-Scholes put */
    double N_minus_d1 = 0.5 * (1.0 + erf(-d1 / sqrt(2.0)));
    double N_minus_d2 = 0.5 * (1.0 + erf(-d2 / sqrt(2.0)));
    double bs_put = strike * exp(-rate*maturity) * N_minus_d2 - spot * exp(-div_yield*maturity) * N_minus_d1;
    
    printf("  Black-Scholes price:    $%.6f\n", bs_put);
    printf("  Absolute error:         $%.6f\n", fabs(put_price - bs_put));
    printf("  Relative error:         %.2f%%\n", 100.0 * fabs(put_price - bs_put) / bs_put);
    printf("\n");
    
    /* Verify put-call parity: C - P = S*exp(-q*T) - K*exp(-r*T) */
    printf("Put-Call Parity Verification:\n");
    printf("--------------------------------------------------\n");
    double parity_lhs = call_price - put_price;
    double parity_rhs = spot * exp(-div_yield * maturity) - strike * exp(-rate * maturity);
    double parity_diff = fabs(parity_lhs - parity_rhs);
    
    printf("  C - P:                  $%.6f\n", parity_lhs);
    printf("  S*exp(-qT) - K*exp(-rT): $%.6f\n", parity_rhs);
    printf("  Difference:             $%.6f\n", parity_diff);
    
    if (parity_diff < 0.01) {
        printf("  Status:                 ✓ PASS (error < $0.01)\n");
    } else {
        printf("  Status:                 ✗ WARNING (error >= $0.01)\n");
    }
    printf("\n");
    
    /* Test at-the-money, in-the-money, out-of-the-money */
    printf("Call Prices at Different Strikes:\n");
    printf("--------------------------------------------------\n");
    double strikes[] = {80.0, 90.0, 100.0, 110.0, 120.0};
    int n_strikes = sizeof(strikes) / sizeof(strikes[0]);
    
    printf("  Strike    Call Price    Intrinsic    Time Value\n");
    printf("  ------    ----------    ---------    ----------\n");
    
    for (int i = 0; i < n_strikes; i++) {
        double K = strikes[i];
        double price = fdp_price_european_call(
            spot, K, rate, div_yield, vol, maturity,
            n_space, n_time
        );
        
        if (price < 0) continue;
        
        double intrinsic = fmax(spot - K * exp(-rate * maturity), 0.0);
        double time_value = price - intrinsic;
        
        printf("  $%-6.0f    $%-9.4f    $%-8.4f    $%-9.4f",
               K, price, intrinsic, time_value);
        
        /* Flag negative time values (impossible!) */
        if (time_value < -0.01) {
            printf("  ← ERROR: Negative time value!");
        }
        printf("\n");
    }
    printf("\n");
    
    /* Test different maturities */
    printf("Call Prices at Different Maturities:\n");
    printf("--------------------------------------------------\n");
    double maturities[] = {0.25, 0.5, 1.0, 2.0, 3.0};
    int n_maturities = sizeof(maturities) / sizeof(maturities[0]);
    
    printf("  Maturity (years)    Call Price      Status\n");
    printf("  ----------------    ----------      ------\n");
    
    for (int i = 0; i < n_maturities; i++) {
        double T = maturities[i];
        double price = fdp_price_european_call(
            spot, strike, rate, div_yield, vol, T,
            n_space, n_time
        );
        
        if (price < 0) {
            printf("  %-18.2f    ERROR\n", T);
            continue;
        }
        
        printf("  %-18.2f    $%-9.6f", T, price);
        
        /* Sanity checks */
        if (price < 0 || price > spot) {
            printf("  ← INVALID");
        } else if (fabs(price) > 1e10) {
            printf("  ← OVERFLOW");
        } else {
            printf("  OK");
        }
        printf("\n");
    }
    printf("\n");
    
    /* Test different volatilities */
    printf("Call Prices at Different Volatilities:\n");
    printf("--------------------------------------------------\n");
    double vols[] = {0.10, 0.15, 0.20, 0.30, 0.40};
    int n_vols = sizeof(vols) / sizeof(vols[0]);
    
    printf("  Volatility    Call Price      Status\n");
    printf("  ----------    ----------      ------\n");
    
    for (int i = 0; i < n_vols; i++) {
        double sigma = vols[i];
        double price = fdp_price_european_call(
            spot, strike, rate, div_yield, sigma, maturity,
            n_space, n_time
        );
        
        if (price < 0) {
            printf("  %-8.0f%%    ERROR\n", sigma * 100.0);
            continue;
        }
        
        printf("  %-8.0f%%    $%-9.6f", sigma * 100.0, price);
        
        /* Sanity checks */
        if (price < 0 || price > spot) {
            printf("  ← INVALID");
        } else if (fabs(price) > 1e10) {
            printf("  ← OVERFLOW");
        } else {
            printf("  OK");
        }
        printf("\n");
    }
    printf("\n");
    
    printf("=================================================\n");
    printf("Example completed successfully!\n");
    printf("=================================================\n");
    
    return 0;
}
