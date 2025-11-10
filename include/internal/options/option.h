/**
 * option.h - Option specification structures
 * 
 * Defines different option types and their payoff functions:
 * - Vanilla (European, American, Bermudan)
 * - Barrier (knock-in, knock-out)
 * - Asian (arithmetic/geometric average)
 * - Lookback (fixed/floating strike)
 */

#ifndef FDPRICING_INTERNAL_OPTION_H
#define FDPRICING_INTERNAL_OPTION_H

#include "fdpricing.h"
#include "internal/core/context.h"

/* Option category (for internal dispatch) */
typedef enum {
    FDP_OPTION_CATEGORY_VANILLA = 0,
    FDP_OPTION_CATEGORY_BARRIER = 1,
    FDP_OPTION_CATEGORY_ASIAN = 2,
    FDP_OPTION_CATEGORY_LOOKBACK = 3,
    FDP_OPTION_CATEGORY_DIGITAL = 4
} fdp_option_category_t;

/* Base option structure */
struct fdp_option_s {
    fdp_context_t* ctx;
    fdp_option_category_t category;
    fdp_option_type_t type;        /* Call or Put */
    fdp_option_style_t style;      /* European, American, Bermudan */
    
    double strike;
    double maturity;
    
    /* Barrier options */
    fdp_barrier_type_t barrier_type;
    double barrier_level;
    double rebate;
    
    /* Bermudan options */
    double* exercise_times;
    int n_exercise_times;
    
    /* Asian options */
    int n_averaging_points;
    
    /* Lookback options */
    int is_fixed_strike;
};

/* ========================================================================
 * Payoff Functions
 * ======================================================================== */

/**
 * Compute terminal payoff for vanilla option
 * Returns max(S - K, 0) for call, max(K - S, 0) for put
 */
double fdp_payoff_vanilla(
    fdp_option_type_t type,
    double S,
    double strike
);

/**
 * Compute terminal payoff for barrier option
 * Must check if barrier was hit during path
 */
double fdp_payoff_barrier(
    fdp_option_type_t type,
    fdp_barrier_type_t barrier_type,
    double S,
    double strike,
    double barrier_level,
    double rebate,
    int barrier_hit
);

/**
 * Compute terminal payoff for digital option
 * Returns 1 if in-the-money, 0 otherwise
 */
double fdp_payoff_digital(
    fdp_option_type_t type,
    double S,
    double strike
);

/* ========================================================================
 * Boundary Conditions
 * ======================================================================== */

/**
 * Get boundary condition at S_min for given option
 */
double fdp_boundary_lower(
    const fdp_option_s* option,
    double t,
    double rate
);

/**
 * Get boundary condition at S_max for given option
 */
double fdp_boundary_upper(
    const fdp_option_s* option,
    double s_max,
    double t,
    double rate
);

/* ========================================================================
 * Early Exercise (American/Bermudan)
 * ======================================================================== */

/**
 * Check if option can be exercised at time t
 */
int fdp_option_can_exercise(
    const fdp_option_s* option,
    double t
);

/**
 * Compute exercise value (intrinsic value) at spot S
 */
double fdp_option_exercise_value(
    const fdp_option_s* option,
    double S
);

#endif /* FDPRICING_INTERNAL_OPTION_H */
