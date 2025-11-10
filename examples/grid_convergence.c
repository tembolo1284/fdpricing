/* examples/grid_convergence.c - Test convergence with grid refinement */

#include "fdpricing.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=================================================\n");
    printf("Grid Convergence Study\n");
    printf("=================================================\n\n");
    
    /* Market parameters */
    double spot = 100.0;
    double strike = 100.0;
    double rate = 0.05;
    double div_yield = 0.0;
    double vol = 0.20;
    double maturity = 1.0;
    
    printf("Testing European call option pricing convergence\n");
    printf("as grid is refined...\n\n");
    
    printf("Market Parameters:\n");
    printf("  Spot:       $%.2f\n", spot);
    printf("  Strike:     $%.2f\n", strike);
    printf("  Rate:       %.2f%%\n", rate * 100.0);
    printf("  Volatility: %.2f%%\n", vol * 100.0);
    printf("  Maturity:   %.2f years\n", maturity);
    printf("\n");
    
    /* Test different grid sizes */
    printf("Grid Size Convergence:\n");
    printf("--------------------------------------------------\n");
    printf("N_space  N_time    Call Price    Change from prev\n");
    printf("-------  ------    ----------    ----------------\n");
    
    int grid_sizes[] = {25, 50, 100, 200, 400};
    int n_sizes = sizeof(grid_sizes) / sizeof(grid_sizes[0]);
    double prev_price = 0.0;
    
    for (int i = 0; i < n_sizes; i++) {
        int n = grid_sizes[i];
        
        double price = fdp_price_european_call(
            spot, strike, rate, div_yield, vol, maturity,
            n, n
        );
        
        if (price < 0) {
            printf("%-7d  %-6d    ERROR\n", n, n);
            continue;
        }
        
        if (i == 0) {
            printf("%-7d  %-6d    $%-9.6f    (baseline)\n", n, n, price);
        } else {
            double change = price - prev_price;
            printf("%-7d  %-6d    $%-9.6f    $%+.6f\n", 
                   n, n, price, change);
        }
        
        prev_price = price;
    }
    printf("\n");
    
    printf("Observation: Price should converge as grid is refined\n");
    printf("             Changes should decrease with finer grids\n");
    printf("\n");
    
    /* Test American put convergence */
    printf("American Put Convergence:\n");
    printf("--------------------------------------------------\n");
    printf("N_space  N_time    Put Price     Change from prev\n");
    printf("-------  ------    ---------     ----------------\n");
    
    prev_price = 0.0;
    
    for (int i = 0; i < n_sizes; i++) {
        int n = grid_sizes[i];
        
        double price = fdp_price_american_put(
            spot, strike, rate, div_yield, vol, maturity,
            n, n
        );
        
        if (price < 0) {
            printf("%-7d  %-6d    ERROR\n", n, n);
            continue;
        }
        
        if (i == 0) {
            printf("%-7d  %-6d    $%-9.6f    (baseline)\n", n, n, price);
        } else {
            double change = price - prev_price;
            printf("%-7d  %-6d    $%-9.6f    $%+.6f\n",
                   n, n, price, change);
        }
        
        prev_price = price;
    }
    printf("\n");
    
    /* Test time step independence */
    printf("Time Step Study (N_space=100):\n");
    printf("--------------------------------------------------\n");
    printf("N_time    Call Price    Change from prev\n");
    printf("------    ----------    ----------------\n");
    
    int time_steps[] = {25, 50, 100, 200, 400};
    int n_time_steps = sizeof(time_steps) / sizeof(time_steps[0]);
    prev_price = 0.0;
    
    for (int i = 0; i < n_time_steps; i++) {
        int nt = time_steps[i];
        
        double price = fdp_price_european_call(
            spot, strike, rate, div_yield, vol, maturity,
            100, nt
        );
        
        if (price < 0) {
            printf("%-6d    ERROR\n", nt);
            continue;
        }
        
        if (i == 0) {
            printf("%-6d    $%-9.6f    (baseline)\n", nt, price);
        } else {
            double change = price - prev_price;
            printf("%-6d    $%-9.6f    $%+.6f\n", nt, price, change);
        }
        
        prev_price = price;
    }
    printf("\n");
    
    printf("=================================================\n");
    printf("Example completed successfully!\n");
    printf("=================================================\n");
    
    return 0;
}
