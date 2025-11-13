/**
 * gbm.c - Geometric Brownian Motion (Black-Scholes) model
 * 
 * SDE: dS = (r - q)S dt + sigma*S dW
 * PDE: dV/dt + (r-q)S*dV/dS + (1/2)*sigma^2*S^2*d^2V/dS^2 - rV = 0
 */

#include "fdpricing.h"
#include "internal/models/model.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"

/* ========================================================================
 * GBM Virtual Table Functions
 * ======================================================================== */

static void gbm_get_coefficients_1d(
    const fdp_model_t* model,
    double S,
    double t,
    double* mu,
    double* sigma,
    double* r)
{
    const fdp_gbm_params_t* params = (const fdp_gbm_params_t*)model->params;
    
    /* For the log-transformed PDE or when discretizing in S-space,
     * we need to return the LOCAL drift and diffusion:
     * 
     * Drift coefficient: μ = (r - q) for ∂V/∂S term
     * Diffusion coefficient: σ = vol for ∂²V/∂S² term
     * 
     * The S and S² factors are handled in the discretization
     */
    
    *mu = (params->rate - params->div_yield);
    *sigma = params->vol;
    
    /* Discount rate */
    *r = params->rate;
    
    (void)S;
    (void)t; 
}

static void gbm_get_coefficients_2d(
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
    /* GBM is 1D only - this should never be called */
    (void)model; (void)S; (void)v; (void)t;
    (void)mu_s; (void)mu_v; (void)sigma_s; (void)sigma_v; (void)rho; (void)r;
}

static void gbm_apply_jump_integral(
    const fdp_model_t* model,
    double* grid_values,
    const double* space_grid,
    int n_space,
    double t)
{
    /* GBM has no jumps - nothing to do */
    (void)model; (void)grid_values; (void)space_grid; (void)n_space; (void)t;
}

static void gbm_destroy(fdp_model_t* model)
{
    if (!model) return;
    
    /* Free parameters */
    if (model->params) {
        fdp_ctx_free(model->ctx, model->params);
    }
}

/* GBM virtual table */
static const fdp_model_vtable_t gbm_vtable = {
    gbm_get_coefficients_1d,
    gbm_get_coefficients_2d,
    gbm_apply_jump_integral,
    gbm_destroy
};

/* ========================================================================
 * GBM Model Creation
 * ======================================================================== */

fdp_model_t* fdp_model_new_gbm_internal(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double vol)
{
    if (!ctx) return NULL;
    
    /* Validate parameters */
    if (vol < 0.0) {
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
    fdp_gbm_params_t* params = FDP_CTX_ALLOC_ARRAY(ctx, fdp_gbm_params_t, 1);
    if (!params) {
        fdp_ctx_free(ctx, model);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Initialize parameters */
    params->rate = rate;
    params->div_yield = div_yield;
    params->vol = vol;
    
    /* Initialize model */
    model->ctx = ctx;
    model->type = FDP_MODEL_GBM;
    model->vtable = &gbm_vtable;
    model->params = params;
    
    return model;
}

/* ========================================================================
 * Public API
 * ======================================================================== */

fdp_model_t* fdp_model_new_gbm(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double vol)
{
    return fdp_model_new_gbm_internal(ctx, rate, div_yield, vol);
}
