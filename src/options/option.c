/**
 * option.c - Option specification and payoff functions
 */

#include "fdpricing.h"
#include "internal/options/option.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"
#include <math.h>

/* ========================================================================
 * Option Creation (Public API)
 * ======================================================================== */

fdp_option_t* fdp_option_new_vanilla(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    fdp_option_style_t style,
    double strike,
    double maturity)
{
    if (!ctx || strike <= 0.0 || maturity <= 0.0) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_VANILLA;
    option->type = type;
    option->style = style;
    option->strike = strike;
    option->maturity = maturity;
    
    /* Initialize unused fields */
    option->barrier_type = FDP_BARRIER_UP_OUT;
    option->barrier_level = 0.0;
    option->rebate = 0.0;
    option->exercise_times = NULL;
    option->n_exercise_times = 0;
    option->n_averaging_points = 0;
    option->is_fixed_strike = 0;
    
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
    if (!ctx || strike <= 0.0 || maturity <= 0.0 || 
        !exercise_times || n_exercise_times <= 0) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Copy exercise times */
    option->exercise_times = FDP_CTX_ALLOC_ARRAY(ctx, double, n_exercise_times);
    if (!option->exercise_times) {
        fdp_ctx_free(ctx, option);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    for (int i = 0; i < n_exercise_times; ++i) {
        option->exercise_times[i] = exercise_times[i];
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_VANILLA;
    option->type = type;
    option->style = FDP_STYLE_BERMUDAN;
    option->strike = strike;
    option->maturity = maturity;
    option->n_exercise_times = n_exercise_times;
    
    /* Initialize unused fields */
    option->barrier_type = FDP_BARRIER_UP_OUT;
    option->barrier_level = 0.0;
    option->rebate = 0.0;
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
    if (!ctx || strike <= 0.0 || maturity <= 0.0 || barrier_level <= 0.0) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
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
    if (!ctx || strike <= 0.0 || maturity <= 0.0 || n_averaging_points <= 0) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
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
    option->barrier_type = FDP_BARRIER_UP_OUT;
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
    if (!ctx || maturity <= 0.0) {
        if (ctx) fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    if (is_fixed_strike && strike <= 0.0) {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    fdp_option_t* option = FDP_CTX_ALLOC_ARRAY(ctx, fdp_option_t, 1);
    if (!option) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    option->ctx = ctx;
    option->category = FDP_OPTION_CATEGORY_LOOKBACK;
    option->type = type;
    option->style = FDP_STYLE_EUROPEAN;
    option->strike = is_fixed_strike ? strike : 0.0;
    option->maturity = maturity;
    option->is_fixed_strike = is_fixed_strike;
    
    /* Initialize unused fields */
    option->barrier_type = FDP_BARRIER_UP_OUT;
    option->barrier_level = 0.0;
    option->rebate = 0.0;
    option->exercise_times = NULL;
    option->n_exercise_times = 0;
    option->n_averaging_points = 0;
    
    return option;
}

void fdp_option_free(fdp_option_t* option)
{
    if (!option) return;
    
    fdp_context_t* ctx = option->ctx;
    
    if (option->exercise_times) {
        fdp_ctx_free(ctx, option->exercise_times);
    }
    
    fdp_ctx_free(ctx, option);
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
        double intrinsic = S - strike;
        return (intrinsic > 0.0) ? intrinsic : 0.0;
    } else {
        double intrinsic = strike - S;
        return (intrinsic > 0.0) ? intrinsic : 0.0;
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
    double vanilla_payoff = fdp_payoff_vanilla(type, S, strike);
    
    switch (barrier_type) {
        case FDP_BARRIER_UP_OUT:
        case FDP_BARRIER_DOWN_OUT:
            /* Knock-out: if barrier hit, get rebate, else vanilla payoff */
            return barrier_hit ? rebate : vanilla_payoff;
            
        case FDP_BARRIER_UP_IN:
        case FDP_BARRIER_DOWN_IN:
            /* Knock-in: if barrier hit, get vanilla payoff, else rebate */
            return barrier_hit ? vanilla_payoff : rebate;
            
        default:
            return 0.0;
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
    const fdp_option_s* option,
    double t,
    double rate)
{
    /* At S = 0 (or very small S) */
    double time_to_maturity = option->maturity - t;
    
    if (option->type == FDP_OPTION_CALL) {
        /* Call is worthless when S = 0 */
        return 0.0;
    } else {
        /* Put approaches discounted strike when S = 0 */
        return option->strike * exp(-rate * time_to_maturity);
    }
}

double fdp_boundary_upper(
    const fdp_option_s* option,
    double s_max,
    double t,
    double rate)
{
    /* At S = S_max (very large S) */
    double time_to_maturity = option->maturity - t;
    
    if (option->type == FDP_OPTION_CALL) {
        /* Call approaches S - K*exp(-r*T) for large S */
        return s_max - option->strike * exp(-rate * time_to_maturity);
    } else {
        /* Put is worthless for large S */
        return 0.0;
    }
}

/* ========================================================================
 * Early Exercise
 * ======================================================================== */

int fdp_option_can_exercise(
    const fdp_option_s* option,
    double t)
{
    if (!option) return 0;
    
    if (option->style == FDP_STYLE_AMERICAN) {
        /* Can exercise any time */
        return 1;
    }
    
    if (option->style == FDP_STYLE_BERMUDAN) {
        /* Check if t matches any exercise time */
        double tolerance = 1e-10;
        for (int i = 0; i < option->n_exercise_times; ++i) {
            if (fabs(t - option->exercise_times[i]) < tolerance) {
                return 1;
            }
        }
        return 0;
    }
    
    /* European: only at maturity */
    return 0;
}

double fdp_option_exercise_value(
    const fdp_option_s* option,
    double S)
{
    if (!option) return 0.0;
    
    /* Exercise value is just the intrinsic value */
    return fdp_payoff_vanilla(option->type, S, option->strike);
}
