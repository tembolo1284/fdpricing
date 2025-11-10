/* examples/american_comparison.c - Compare American and European options */

#include "fdpricing.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=================================================\n");
    printf("American vs European Option Comparison\n");
    printf("=================================================\n\n");
    
    /* Market parameters */
    double spot = 100.0;
    double strike = 100.0;
    double rate = 0.05;
    double div_yield = 0.0;
    double vol = 0.20;
    double maturity = 1.0;
    
    /* Grid parameters */
    int n_space = 150;  /* Use finer grid for American options */
    int n_time = 150;
    
    printf("Market Parameters:\n");
    printf("  Spot price:             $%.2f\n", spot);
    printf("  Strike price:           $%.2f\n", strike);
    printf("  Risk-free rate:         %.2f%%\n", rate * 100.0);
    printf("  Dividend yield:         %.2f%%\n", div_yield * 100.0);
    printf("  Volatility:             %.2f%%\n", vol * 100.0);
    printf("  Time to maturity:       %.2f years\n", maturity);
    printf("\n");
    
    /* Price American and European calls */
    printf("CALL OPTIONS:\n");
    printf("--------------------------------------------------\n");
    
    double euro_call = fdp_price_european_call(
        spot, strike, rate, div_yield, vol, maturity,
        n_space, n_time
    );
    
    double amer_call = fdp_price_american_call(
        spot, strike, rate, div_yield, vol, maturity,
        n_space, n_time
    );
    
    if (euro_call < 0 || amer_call < 0) {
        printf("  ERROR: Call pricing failed!\n");
        return 1;
    }
    
    printf("  European call:          $%.6f\n", euro_call);
    printf("  American call:          $%.6f\n", amer_call);
    printf("  Early exercise premium: $%.6f (%.2f%%)\n",
           amer_call - euro_call,
           100.0 * (amer_call - euro_call) / euro_call);
    printf("\n");
    
    printf("Note: For calls with no dividends, American = European\n");
    printf("      (no early exercise benefit)\n");
    printf("\n");
    
    /* Price American and European puts */
    printf("PUT OPTIONS:\n");
    printf("--------------------------------------------------\n");
    
    double euro_put = fdp_price_european_put(
        spot, strike, rate, div_yield, vol, maturity,
        n_space, n_time
    );
    
    double amer_put = fdp_price_american_put(
        spot, strike, rate, div_yield, vol, maturity,
        n_space, n_time
    );
    
    if (euro_put < 0 || amer_put < 0) {
        printf("  ERROR: Put pricing failed!\n");
        return 1;
    }
    
    printf("  European put:           $%.6f\n", euro_put);
    printf("  American put:           $%.6f\n", amer_put);
    printf("  Early exercise premium: $%.6f (%.2f%%)\n",
           amer_put - euro_put,
           100.0 * (amer_put - euro_put) / euro_put);
    printf("\n");
    
    printf("Note: American puts have early exercise value\n");
    printf("      when deep in-the-money\n");
    printf("\n");
    
    /* Test early exercise premium at different strikes */
    printf("Early Exercise Premium vs Strike (Puts):\n");
    printf("--------------------------------------------------\n");
    printf("Strike    European    American    Premium    Premium%%\n");
    printf("------    --------    --------    -------    --------\n");
    
    double strikes[] = {80.0, 90.0, 100.0, 110.0, 120.0, 130.0};
    int n_strikes = sizeof(strikes) / sizeof(strikes[0]);
    
    for (int i = 0; i < n_strikes; i++) {
        double K = strikes[i];
        
        double ep = fdp_price_european_put(
            spot, K, rate, div_yield, vol, maturity,
            n_space, n_time
        );
        
        double ap = fdp_price_american_put(
            spot, K, rate, div_yield, vol, maturity,
            n_space, n_time
        );
        
        if (ep < 0 || ap < 0) continue;
        
        double premium = ap - ep;
        double premium_pct = (ep > 0) ? 100.0 * premium / ep : 0.0;
        
        printf("$%-5.0f    $%-7.4f    $%-7.4f    $%-6.4f    %6.2f%%\n",
               K, ep, ap, premium, premium_pct);
    }
    printf("\n");
    
    printf("Observation: Early exercise premium increases for\n");
    printf("             deep in-the-money puts (high strike)\n");
    printf("\n");
    
    /* Test early exercise premium at different maturities */
    printf("Early Exercise Premium vs Maturity (Puts, K=110):\n");
    printf("--------------------------------------------------\n");
    printf("Maturity    European    American    Premium    Premium%%\n");
    printf("--------    --------    --------    -------    --------\n");
    
    double K_test = 110.0;
    double maturities[] = {0.25, 0.5, 1.0, 2.0, 3.0};
    int n_maturities = sizeof(maturities) / sizeof(maturities[0]);
    
    for (int i = 0; i < n_maturities; i++) {
        double T = maturities[i];
        
        double ep = fdp_price_european_put(
            spot, K_test, rate, div_yield, vol, T,
            n_space, n_time
        );
        
        double ap = fdp_price_american_put(
            spot, K_test, rate, div_yield, vol, T,
            n_space, n_time
        );
        
        if (ep < 0 || ap < 0) continue;
        
        double premium = ap - ep;
        double premium_pct = (ep > 0) ? 100.0 * premium / ep : 0.0;
        
        printf("%-6.2f yr    $%-7.4f    $%-7.4f    $%-6.4f    %6.2f%%\n",
               T, ep, ap, premium, premium_pct);
    }
    printf("\n");
    
    /* Test with dividends - American calls become valuable */
    printf("EFFECT OF DIVIDENDS:\n");
    printf("--------------------------------------------------\n");
    
    double div_test = 0.03;  /* 3% dividend yield */
    
    printf("With dividend yield = %.1f%%:\n\n", div_test * 100.0);
    
    double euro_call_div = fdp_price_european_call(
        spot, strike, rate, div_test, vol, maturity,
        n_space, n_time
    );
    
    double amer_call_div = fdp_price_american_call(
        spot, strike, rate, div_test, vol, maturity,
        n_space, n_time
    );
    
    if (euro_call_div >= 0 && amer_call_div >= 0) {
        printf("  European call:          $%.6f\n", euro_call_div);
        printf("  American call:          $%.6f\n", amer_call_div);
        printf("  Early exercise premium: $%.6f (%.2f%%)\n",
               amer_call_div - euro_call_div,
               100.0 * (amer_call_div - euro_call_div) / euro_call_div);
        printf("\n");
        printf("Note: With dividends, American calls may be exercised\n");
        printf("      early to capture dividend payments\n");
    }
    printf("\n");
    
    printf("=================================================\n");
    printf("Example completed successfully!\n");
    printf("=================================================\n");
    
    return 0;
}
