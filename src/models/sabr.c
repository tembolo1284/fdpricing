/**
 * sabr.c - SABR stochastic volatility model
 * 
 * SABR (Stochastic Alpha Beta Rho) model widely used for interest rates and FX.
 * 
 * System of SDEs:
 *   dF = alpha * F^beta * dW1
 *   dalpha = nu * alpha * dW2
 *   Cor(dW1, dW2) = rho
 * 
 * Where:
 *   F = forward price (or forward rate)
 *   alpha = stochastic volatility parameter
 *   beta = CEV exponent (0 = normal vol, 1 = lognormal vol)
 *   nu = vol-of-vol
 *   rho = correlation between forward and volatility
 * 
 * 2D PDE for V(F, alpha, t):
 *   dV/dt + (1/2)*alpha^2*F^(2*beta)*d^2V/dF^2 
 *   + (1/2)*nu^2*alpha^2*d^2V/dalpha^2
 *   + rho*nu*alpha^2*F^beta*d^2V/dFdalpha - r*V = 0
 * 
 * Note: Forward has no drift under T-forward measure (martingale)
 */

#include "fdpricing.h"
#include "internal/models/model.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"
#include <math.h>

/* ========================================================================
 * SABR Virtual Table Functions
 * ======================================================================== */

static void sabr_get_coefficients_1d(
    const fdp_model_t* model,
    double S,
    double t,
    double* mu,
    double* sigma,
    double* r)
{
    /* SABR is a 2D model - for 1D approximation use alpha at current level */
    const fdp_sabr_params_t* params = (const fdp_sabr_params_t*)model->params;
    
    /* No drift for forward (martingale) */
    *mu = 0.0;
    
    /* Effective volatility: sigma_eff = alpha * F^beta */
    double F_beta = pow(S, params->beta);
    *sigma = params->alpha * F_beta * S;  /* Need to multiply by S for PDE form */
    
    /* Discount rate */
    *r = params->rate;
    
    (void)t;
}

static void sabr_get_coefficients_2d(
    const fdp_model_t* model,
    double F,
    double alpha,
    double t,
    double* mu_f,
    double* mu_alpha,
    double* sigma_f,
    double* sigma_alpha,
    double* rho,
    double* r)
{
    const fdp_sabr_params_t* params = (const fdp_sabr_params_t*)model->params;
    
    /* Ensure positive values */
    double F_pos = (F > 0.0) ? F : 1e-10;
    double alpha_pos = (alpha > 0.0) ? alpha : 1e-10;
    
    /* Compute F^beta for diffusion term */
    double F_beta = pow(F_pos, params->beta);
    
    /* Drift in F: zero (forward is martingale under T-forward measure) */
    *mu_f = 0.0;
    
    /* Drift in alpha: zero (alpha is driftless martingale) */
    *mu_alpha = 0.0;
    
    /* Diffusion in F: sigma_f = alpha * F^beta */
    *sigma_f = alpha_pos * F_beta * F_pos;  /* Multiply by F for PDE coefficient form */
    
    /* Diffusion in alpha: sigma_alpha = nu * alpha */
    *sigma_alpha = params->nu * alpha_pos;
    
    /* Correlation */
    *rho = params->rho;
    
    /* Discount rate */
    *r = params->rate;
    
    (void)t;  /* Time-homogeneous */
}

static void sabr_apply_jump_integral(
    const fdp_model_t* model,
    double* grid_values,
    const double* space_grid,
    int n_space,
    double t)
{
    /* SABR has no jumps */
    (void)model; (void)grid_values; (void)space_grid; (void)n_space; (void)t;
}

static void sabr_destroy(fdp_model_t* model)
{
    if (!model) return;
    
    if (model->params) {
        fdp_ctx_free(model->ctx, model->params);
    }
}

/* SABR virtual table */
static const fdp_model_vtable_t sabr_vtable = {
    sabr_get_coefficients_1d,
    sabr_get_coefficients_2d,
    sabr_apply_jump_integral,
    sabr_destroy
};

/* ========================================================================
 * SABR Model Creation
 * ======================================================================== */

fdp_model_t* fdp_model_new_sabr_internal(
    fdp_context_t* ctx,
    double rate,
    double forward,
    double alpha,
    double beta,
    double nu,
    double rho)
{
    if (!ctx) return NULL;
    
    /* Validate parameters */
    if (forward <= 0.0 || alpha <= 0.0 || nu < 0.0) {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    /* Beta should be in [0, 1] typically, but we'll allow any value */
    /* Common values: 0 (normal), 0.5 (CIR-like), 1 (lognormal) */
    
    if (rho < -1.0 || rho > 1.0) {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    /* Allocate model structure */
    fdp_model_t* model = FDP_CTX_ALLOC_ARRAY(ctx, fdp_model_t, 1);
    if (!model) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Allocate parameters */
    fdp_sabr_params_t* params = FDP_CTX_ALLOC_ARRAY(ctx, fdp_sabr_params_t, 1);
    if (!params) {
        fdp_ctx_free(ctx, model);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Initialize parameters */
    params->rate = rate;
    params->forward = forward;
    params->alpha = alpha;
    params->beta = beta;
    params->nu = nu;
    params->rho = rho;
    
    /* Initialize model */
    model->ctx = ctx;
    model->type = FDP_MODEL_SABR;
    model->vtable = &sabr_vtable;
    model->params = params;
    
    return model;
}

/* ========================================================================
 * Public API
 * ======================================================================== */

fdp_model_t* fdp_model_new_sabr(
    fdp_context_t* ctx,
    double rate,
    double forward,
    double alpha,
    double beta,
    double nu,
    double rho)
{
    return fdp_model_new_sabr_internal(ctx, rate, forward, alpha, beta, nu, rho);
}
