/**
 * merton_jump.c - Merton jump-diffusion model
 * 
 * SDE: dS = (r - q - lambda*k)S dt + sigma*S dW + S(e^J - 1)dN
 * where N is Poisson process with intensity lambda
 *       J ~ Normal(mu_j, sigma_j^2)
 *       k = E[e^J - 1] = exp(mu_j + sigma_j^2/2) - 1
 * 
 * PIDE: dV/dt + (r-q-lambda*k)S*dV/dS + (1/2)*sigma^2*S^2*d^2V/dS^2 - rV
 *       + lambda * Integral[ V(S*e^y) * phi(y) dy ] = 0
 * 
 * where phi(y) is the jump size density (lognormal)
 */

#include "fdpricing.h"
#include "internal/models/model.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========================================================================
 * Merton Virtual Table Functions
 * ======================================================================== */

static void merton_get_coefficients_1d(
    const fdp_model_t* model,
    double S,
    double t,
    double* mu,
    double* sigma,
    double* r)
{
    const fdp_merton_params_t* params = (const fdp_merton_params_t*)model->params;
    
    /* Expected jump size: k = E[e^J - 1] = exp(mu_j + sigma_j^2/2) - 1 */
    double k = exp(params->mu_j + 0.5 * params->sigma_j * params->sigma_j) - 1.0;
    
    /* Drift compensated for jumps: mu = (r - q - lambda*k)S */
    *mu = (params->rate - params->div_yield - params->lambda * k) * S;
    
    /* Diffusion: sigma = vol * S */
    *sigma = params->vol * S;
    
    /* Discount rate */
    *r = params->rate;
    
    (void)t;
}

static void merton_get_coefficients_2d(
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
    /* Merton is 1D only */
    (void)model; (void)S; (void)v; (void)t;
    (void)mu_s; (void)mu_v; (void)sigma_s; (void)sigma_v; (void)rho; (void)r;
}

static void merton_apply_jump_integral(
    const fdp_model_t* model,
    double* grid_values,
    const double* space_grid,
    int n_space,
    double t)
{
    const fdp_merton_params_t* params = (const fdp_merton_params_t*)model->params;
    
    (void)t;  /* Time-homogeneous */
    
    if (params->lambda == 0.0) return;  /* No jumps */
    
    /* We need to compute: lambda * Integral[ V(S*e^y) * phi(y) dy ]
     * where phi(y) = (1/sqrt(2*pi*sigma_j^2)) * exp(-(y-mu_j)^2/(2*sigma_j^2))
     * 
     * We'll approximate this with numerical quadrature.
     * Use a truncated integral: integrate from mu_j - 4*sigma_j to mu_j + 4*sigma_j
     */
    
    const int n_quad = 20;  /* Number of quadrature points */
    double y_min = params->mu_j - 4.0 * params->sigma_j;
    double y_max = params->mu_j + 4.0 * params->sigma_j;
    double dy = (y_max - y_min) / (double)n_quad;
    
    /* Allocate temporary array for jump contributions */
    double* jump_contrib = FDP_CTX_CALLOC_ARRAY(model->ctx, double, n_space);
    if (!jump_contrib) return;
    
    /* Compute jump density normalization */
    double sqrt_2pi_sigma = sqrt(2.0 * M_PI) * params->sigma_j;
    double two_sigma_sq = 2.0 * params->sigma_j * params->sigma_j;
    
    /* For each grid point, integrate over jump sizes */
    for (int i = 0; i < n_space; ++i) {
        double S = space_grid[i];
        double integral = 0.0;
        
        /* Numerical integration using midpoint rule */
        for (int j = 0; j < n_quad; ++j) {
            double y = y_min + ((double)j + 0.5) * dy;
            
            /* Jump density: phi(y) */
            double deviation = y - params->mu_j;
            double phi = exp(-deviation * deviation / two_sigma_sq) / sqrt_2pi_sigma;
            
            /* New stock price after jump: S_new = S * e^y */
            double S_new = S * exp(y);
            
            /* Interpolate V(S_new) from grid */
            double V_new = 0.0;
            
            if (S_new <= space_grid[0]) {
                V_new = grid_values[0];
            } else if (S_new >= space_grid[n_space - 1]) {
                V_new = grid_values[n_space - 1];
            } else {
                /* Linear interpolation */
                int k = 0;
                while (k < n_space - 1 && space_grid[k + 1] < S_new) {
                    ++k;
                }
                
                double S0 = space_grid[k];
                double S1 = space_grid[k + 1];
                double V0 = grid_values[k];
                double V1 = grid_values[k + 1];
                
                double alpha = (S_new - S0) / (S1 - S0);
                V_new = V0 + alpha * (V1 - V0);
            }
            
            integral += V_new * phi * dy;
        }
        
        jump_contrib[i] = params->lambda * integral;
    }
    
    /* Add jump contributions to grid values */
    for (int i = 0; i < n_space; ++i) {
        grid_values[i] += jump_contrib[i];
    }
    
    fdp_ctx_free(model->ctx, jump_contrib);
}

static void merton_destroy(fdp_model_t* model)
{
    if (!model) return;
    
    if (model->params) {
        fdp_ctx_free(model->ctx, model->params);
    }
}

/* Merton virtual table */
static const fdp_model_vtable_t merton_vtable = {
    merton_get_coefficients_1d,
    merton_get_coefficients_2d,
    merton_apply_jump_integral,
    merton_destroy
};

/* ========================================================================
 * Merton Model Creation
 * ======================================================================== */

fdp_model_t* fdp_model_new_merton_internal(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double vol,
    double lambda,
    double mu_j,
    double sigma_j)
{
    if (!ctx) return NULL;
    
    /* Validate parameters */
    if (vol < 0.0 || lambda < 0.0 || sigma_j < 0.0) {
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
    fdp_merton_params_t* params = FDP_CTX_ALLOC_ARRAY(ctx, fdp_merton_params_t, 1);
    if (!params) {
        fdp_ctx_free(ctx, model);
        fdp_ctx_set_error(ctx, FDP_ERROR_ALLOCATION);
        return NULL;
    }
    
    /* Initialize parameters */
    params->rate = rate;
    params->div_yield = div_yield;
    params->vol = vol;
    params->lambda = lambda;
    params->mu_j = mu_j;
    params->sigma_j = sigma_j;
    
    /* Initialize model */
    model->ctx = ctx;
    model->type = FDP_MODEL_MERTON_JUMP;
    model->vtable = &merton_vtable;
    model->params = params;
    
    return model;
}

/* ========================================================================
 * Public API
 * ======================================================================== */

fdp_model_t* fdp_model_new_merton_jump(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double vol,
    double lambda,
    double mu_j,
    double sigma_j)
{
    return fdp_model_new_merton_internal(ctx, rate, div_yield, vol,
                                         lambda, mu_j, sigma_j);
}
