/**
 * heston.c - Heston stochastic volatility model
 * 
 * System of SDEs:
 *   dS = (r - q)S dt + sqrt(v)*S dW1
 *   dv = kappa*(theta - v) dt + sigma*sqrt(v) dW2
 *   Cor(dW1, dW2) = rho
 * 
 * 2D PDE for V(S, v, t):
 *   dV/dt + (r-q)S*dV/dS + kappa*(theta-v)*dV/dv 
 *   + (1/2)*v*S^2*d^2V/dS^2 + (1/2)*sigma^2*v*d^2V/dv^2
 *   + rho*sigma*v*S*d^2V/dSdv - rV = 0
 */

#include "fdpricing.h"
#include "internal/models/model.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"
#include <math.h>

/* ========================================================================
 * Heston Virtual Table Functions
 * ======================================================================== */

static void heston_get_coefficients_1d(
    const fdp_model_t* model,
    double S,
    double t,
    double* mu,
    double* sigma,
    double* r)
{
    /* Heston is a 2D model - this should never be called for full pricing */
    /* However, we can provide effective 1D coefficients at v0 for testing */
    const fdp_heston_params_t* params = (const fdp_heston_params_t*)model->params;
    
    *mu = (params->rate - params->div_yield) * S;
    *sigma = sqrt(params->v0) * S;  /* Use initial variance */
    *r = params->rate;
    
    (void)t;
}

static void heston_get_coefficients_2d(
    const fdp_model_t* model,
    double S,
    double v,
    double t,
    double* mu_s,
    double* mu_v,
    double* sigma_s,
    double* sigma_v,
    double* rho,
    double* r)
{
    const fdp_heston_params_t* params = (const fdp_heston_params_t*)model->params;
    
    /* Ensure variance is non-negative (Feller condition handling) */
    double v_pos = (v > 0.0) ? v : 0.0;
    double sqrt_v = sqrt(v_pos);
    
    /* Drift in S: mu_s = (r - q)S */
    *mu_s = (params->rate - params->div_yield) * S;
    
    /* Drift in v: mu_v = kappa*(theta - v) */
    *mu_v = params->kappa * (params->theta - v);
    
    /* Diffusion in S: sigma_s = sqrt(v)*S */
    *sigma_s = sqrt_v * S;
    
    /* Diffusion in v: sigma_v = sigma*sqrt(v) */
    *sigma_v = params->sigma * sqrt_v;
    
    /* Correlation */
    *rho = params->rho;
    
    /* Discount rate */
    *r = params->rate;
    
    (void)t;  /* Time-homogeneous */
}

static void heston_apply_jump_integral(
    const fdp_model_t* model,
    double* grid_values,
    const double* space_grid,
    int n_space,
    double t)
{
    /* Heston has no jumps */
    (void)model; (void)grid_values; (void)space_grid; (void)n_space; (void)t;
}

static void heston_destroy(fdp_model_t* model)
{
    if (!model) return;
    
    if (model->params) {
        fdp_ctx_free(model->ctx, model->params);
    }
}

/* Heston virtual table */
static const fdp_model_vtable_t heston_vtable = {
    heston_get_coefficients_1d,
    heston_get_coefficients_2d,
    heston_apply_jump_integral,
    heston_destroy
};

/* ========================================================================
 * Heston Model Creation
 * ======================================================================== */

fdp_model_t* fdp_model_new_heston_internal(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double kappa,
    double theta,
    double sigma,
    double rho,
    double v0)
{
    if (!ctx) return NULL;
    
    /* Validate parameters */
    if (kappa < 0.0 || theta < 0.0 || sigma < 0.0 || v0 < 0.0) {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    if (rho < -1.0 || rho > 1.0) {
        fdp_ctx_set_error(ctx, FDP_ERROR_INVALID_PARAM);
        return NULL;
    }
    
    /* Check Feller condition: 2*kappa*theta > sigma^2 */
    /* If violated, variance can reach zero - still solvable but needs care */
    
    /* Allocate model structure */
    fdp_model_t* model = FDP_CTX_ALLOC_ARRAY(ctx, fdp_model_t, 1);
    if (!model) {
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Allocate parameters */
    fdp_heston_params_t* params = FDP_CTX_ALLOC_ARRAY(ctx, fdp_heston_params_t, 1);
    if (!params) {
        fdp_ctx_free(ctx, model);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Initialize parameters */
    params->rate = rate;
    params->div_yield = div_yield;
    params->kappa = kappa;
    params->theta = theta;
    params->sigma = sigma;
    params->rho = rho;
    params->v0 = v0;
    
    /* Initialize model */
    model->ctx = ctx;
    model->type = FDP_MODEL_HESTON;
    model->vtable = &heston_vtable;
    model->params = params;
    
    return model;
}

/* ========================================================================
 * Public API
 * ======================================================================== */

fdp_model_t* fdp_model_new_heston(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double kappa,
    double theta,
    double sigma,
    double rho,
    double v0)
{
    return fdp_model_new_heston_internal(ctx, rate, div_yield, kappa, 
                                         theta, sigma, rho, v0);
}
