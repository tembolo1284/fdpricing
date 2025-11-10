/* src/options/option.c - Option utilities and payoff functions */

#include "fdpricing.h"
#include "internal/core/context.h"
#include "internal/options/option.h"
#include "internal/utils/allocator.h"
#include <math.h>
#include <stdlib.h>

/* ========================================================================
 * Option Creation and Destruction
 * ======================================================================== */

fdp_option_t* fdp_option_new_vanilla(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    fdp_option_style_t style,
    double strike,
    double maturity)
{
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_VANILLA;
    option->type = type;
    option->style = style;
    option->strike = strike;
    option->maturity = maturity;
    
    /* Initialize unused fields */
    option->barrier_type = FDP_BARRIER_NONE;
    option->barrier_level = 0.0;
    option->rebate = 0.0;
    option->exercise_times = NULL;
    option->n_exercise_times = 0;
    option->n_averaging_points = 0;
    option->is_fixed_strike = 0;
    
    return option;
}

fdp_option_t* fdp_option_new_barrier(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,
    double maturity,
    fdp_barrier_type_t barrier_type,
    double barrier_level,
    double rebate)
{
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_BARRIER;
    option->type = type;
    option->style = FDP_STYLE_EUROPEAN;  /* Barriers are typically European */
    option->strike = strike;
    option->maturity = maturity;
    option->barrier_type = barrier_type;
    option->barrier_level = barrier_level;
    option->rebate = rebate;
    
    /* Initialize unused fields */
    option->exercise_times = NULL;
    option->n_exercise_times = 0;
    option->n_averaging_points = 0;
    option->is_fixed_strike = 0;
    
    return option;
}

fdp_option_t* fdp_option_new_asian(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,
    double maturity,
    int n_averaging_points)
{
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_ASIAN;
    option->type = type;
    option->style = FDP_STYLE_EUROPEAN;
    option->strike = strike;
    option->maturity = maturity;
    option->n_averaging_points = n_averaging_points;
    
    /* Initialize unused fields */
    option->barrier_type = FDP_BARRIER_NONE;
    option->barrier_level = 0.0;
    option->rebate = 0.0;
    option->exercise_times = NULL;
    option->n_exercise_times = 0;
    option->is_fixed_strike = 0;
    
    return option;
}

fdp_option_t* fdp_option_new_lookback(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,
    double maturity,
    int is_fixed_strike)
{
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_LOOKBACK;
    option->type = type;
    option->style = FDP_STYLE_EUROPEAN;
    option->strike = strike;
    option->maturity = maturity;
    option->is_fixed_strike = is_fixed_strike;
    
    /* Initialize unused fields */
    option->barrier_type = FDP_BARRIER_NONE;
    option->barrier_level = 0.0;
    option->rebate = 0.0;
    option->exercise_times = NULL;
    option->n_exercise_times = 0;
    option->n_averaging_points = 0;
    
    return option;
}

fdp_option_t* fdp_option_new_bermudan(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,
    double maturity,
    const double* exercise_times,
    int n_exercise_times)
{
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Copy exercise times */
    double* times = FDP_CTX_ALLOC_ARRAY(ctx, double, n_exercise_times);
    if (!times) {
        FDP_CTX_FREE(ctx, option);
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    for (int i = 0; i < n_exercise_times; ++i) {
        times[i] = exercise_times[i];
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_VANILLA;
    option->type = type;
    option->style = FDP_STYLE_BERMUDAN;
    option->strike = strike;
    option->maturity = maturity;
    option->exercise_times = times;
    option->n_exercise_times = n_exercise_times;
    
    /* Initialize unused fields */
    option->barrier_type = FDP_BARRIER_NONE;
    option->barrier_level = 0.0;
    option->rebate = 0.0;
    option->n_averaging_points = 0;
    option->is_fixed_strike = 0;
    
    return option;
}

void fdp_option_free(fdp_option_t* option)
{
    if (!option) return;
    
    /* Free exercise times if allocated */
    if (option->exercise_times) {
        FDP_CTX_FREE(option->ctx, option->exercise_times);
    }
    
    /* Free the option itself */
    FDP_CTX_FREE(option->ctx, option);
}

/* ========================================================================
 * Payoff Functions
 * ======================================================================== */

double fdp_payoff_vanilla(
    fdp_option_type_t type,
    double S,
    double strike)
{
    if (type == FDP_OPTION_CALL) {
        return fmax(S - strike, 0.0);
    } else {
        return fmax(strike - S, 0.0);
    }
}

double fdp_payoff_barrier(
    fdp_option_type_t type,
    fdp_barrier_type_t barrier_type,
    double S,
    double strike,
    double barrier_level,
    double rebate,
    int barrier_hit)
{
    (void)barrier_level;
    double vanilla_payoff = fdp_payoff_vanilla(type, S, strike);
    
    switch (barrier_type) {
        case FDP_BARRIER_UP_OUT:
        case FDP_BARRIER_DOWN_OUT:
            /* Knock-out: pay vanilla if barrier NOT hit, else rebate */
            return barrier_hit ? rebate : vanilla_payoff;
            
        case FDP_BARRIER_UP_IN:
        case FDP_BARRIER_DOWN_IN:
            /* Knock-in: pay vanilla if barrier WAS hit, else rebate */
            return barrier_hit ? vanilla_payoff : rebate;
            
        default:
            return vanilla_payoff;
    }
}

double fdp_payoff_digital(
    fdp_option_type_t type,
    double S,
    double strike)
{
    if (type == FDP_OPTION_CALL) {
        return (S > strike) ? 1.0 : 0.0;
    } else {
        return (S < strike) ? 1.0 : 0.0;
    }
}

/* ========================================================================
 * Boundary Conditions
 * ======================================================================== */

double fdp_boundary_lower(
    const fdp_option_t* option,  /* FIXED: was fdp_option_s* */
    double t,
    double rate)
{
    /* At S = 0, typically option values are known */
    double time_to_maturity = option->maturity - t;
    
    if (option->type == FDP_OPTION_CALL) {
        /* Call worth 0 at S=0 */
        return 0.0;
    } else {
        /* Put worth K*exp(-r*(T-t)) at S=0 */
        return option->strike * exp(-rate * time_to_maturity);
    }
}

double fdp_boundary_upper(
    const fdp_option_t* option,  /* FIXED: was fdp_option_s* */
    double s_max,
    double t,
    double rate)
{
    /* At S = S_max (very large), option behaves like stock or strike */
    double time_to_maturity = option->maturity - t;
    
    if (option->type == FDP_OPTION_CALL) {
        /* Call worth S - K*exp(-r*(T-t)) */
        return s_max - option->strike * exp(-rate * time_to_maturity);
    } else {
        /* Put worth 0 at very large S */
        return 0.0;
    }
}

/* ========================================================================
 * Early Exercise (American/Bermudan)
 * ======================================================================== */

int fdp_option_can_exercise(
    const fdp_option_t* option,  /* FIXED: was fdp_option_s* */
    double t)
{
    /* American options can always exercise */
    if (option->style == FDP_STYLE_AMERICAN) {
        return 1;
    }
    
    /* Bermudan options can only exercise at specific times */
    if (option->style == FDP_STYLE_BERMUDAN) {
        double tolerance = 1e-8;
        /* Check if current time matches any exercise date */
        for (int i = 0; i < option->n_exercise_times; ++i) {
            if (fabs(t - option->exercise_times[i]) < tolerance) {
                return 1;
            }
        }
    }
    
    /* European options cannot exercise early */
    return 0;
}

double fdp_option_exercise_value(
    const fdp_option_t* option,  /* FIXED: was fdp_option_s* */
    double S)
{
    /* Exercise value is just the intrinsic value */
    /* For vanilla options: max(S-K, 0) for calls, max(K-S, 0) for puts */
    return fdp_payoff_vanilla(option->type, S, option->strike);
}
