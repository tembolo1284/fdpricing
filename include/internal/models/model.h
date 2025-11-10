/**
 * model.h - Model interface and structures
 * 
 * Defines the interface for different stochastic models:
 * - GBM (Black-Scholes): dS = rS dt + sigma*S dW
 * - Heston: dS = rS dt + sqrt(v)*S dW1, dv = kappa*(theta-v)dt + sigma*sqrt(v) dW2
 * - Jump models: Add jump integral terms
 */

#ifndef FDPRICING_INTERNAL_MODEL_H
#define FDPRICING_INTERNAL_MODEL_H

#include "fdpricing.h"
#include "internal/core/context.h"

/* Forward declarations - don't redefine typedef, it's in fdpricing.h */
struct fdp_model_s;

/* Model virtual table for polymorphism */
typedef struct {
    /**
     * Get PDE coefficients at point (S, t) for 1D models
     * PDE form: dV/dt + mu(S,t)*dV/dS + (1/2)*sigma^2(S,t)*d^2V/dS^2 - r(S,t)*V = 0
     */
    void (*get_coefficients_1d)(
        const fdp_model_t* model,
        double S,
        double t,
        double* mu,      /* Drift coefficient */
        double* sigma,   /* Diffusion coefficient */
        double* r        /* Discount rate */
    );
    
    /**
     * Get PDE coefficients at point (S, v, t) for 2D models (stochastic vol)
     * PDE form: dV/dt + mu_s*dV/dS + mu_v*dV/dv + (1/2)*sigma_s^2*d^2V/dS^2 
     *           + (1/2)*sigma_v^2*d^2V/dv^2 + rho*sigma_s*sigma_v*d^2V/dSdv - r*V = 0
     */
    void (*get_coefficients_2d)(
        const fdp_model_t* model,
        double S,
        double v,
        double t,
        double* mu_s,    /* Drift in S */
        double* mu_v,    /* Drift in v */
        double* sigma_s, /* Diffusion in S */
        double* sigma_v, /* Diffusion in v */
        double* rho,     /* Correlation */
        double* r        /* Discount rate */
    );
    
    /**
     * Apply jump integral (for jump-diffusion models)
     * Adds integral(V(S*e^y) - V(S)) * k(y) dy term to the PDE
     * 
     * Modifies grid_values in place by adding jump contributions
     */
    void (*apply_jump_integral)(
        const fdp_model_t* model,
        double* grid_values,
        const double* space_grid,
        int n_space,
        double t
    );
    
    /**
     * Cleanup model-specific data
     */
    void (*destroy)(fdp_model_t* model);
    
} fdp_model_vtable_t;

/* Base model structure */
struct fdp_model_s {
    fdp_context_t* ctx;
    fdp_model_type_t type;
    const fdp_model_vtable_t* vtable;
    void* params;  /* Model-specific parameters */
};

/* ========================================================================
 * GBM (Black-Scholes) Model
 * ======================================================================== */

typedef struct {
    double rate;       /* Risk-free rate r */
    double div_yield;  /* Dividend yield q */
    double vol;        /* Volatility sigma */
} fdp_gbm_params_t;

fdp_model_t* fdp_model_new_gbm_internal(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double vol
);

/* ========================================================================
 * Heston Stochastic Volatility Model
 * ======================================================================== */

typedef struct {
    double rate;       /* Risk-free rate r */
    double div_yield;  /* Dividend yield q */
    double kappa;      /* Mean reversion speed */
    double theta;      /* Long-term variance */
    double sigma;      /* Vol of vol (xi) */
    double rho;        /* Correlation rho */
    double v0;         /* Initial variance */
} fdp_heston_params_t;

fdp_model_t* fdp_model_new_heston_internal(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double kappa,
    double theta,
    double sigma,
    double rho,
    double v0
);

/* ========================================================================
 * Merton Jump-Diffusion Model
 * ======================================================================== */

typedef struct {
    double rate;       /* Risk-free rate r */
    double div_yield;  /* Dividend yield q */
    double vol;        /* Diffusion volatility sigma */
    double lambda;     /* Jump intensity lambda */
    double mu_j;       /* Jump mean (in log) */
    double sigma_j;    /* Jump std dev (in log) */
} fdp_merton_params_t;

fdp_model_t* fdp_model_new_merton_internal(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double vol,
    double lambda,
    double mu_j,
    double sigma_j
);

/* ========================================================================
 * Model Utilities
 * ======================================================================== */

/**
 * Check if model is 1D or 2D
 */
int fdp_model_is_1d(const fdp_model_t* model);
int fdp_model_is_2d(const fdp_model_t* model);

/**
 * Check if model has jumps
 */
int fdp_model_has_jumps(const fdp_model_t* model);

/* ========================================================================
 * SABR Stochastic Volatility Model (for forwards/rates)
 * ======================================================================== */

typedef struct {
    double rate;       /* Risk-free rate r (for discounting) */
    double forward;    /* Forward price F */
    double alpha;      /* Initial stochastic volatility parameter */
    double beta;       /* CEV exponent (0 = normal, 0.5 = CIR, 1 = lognormal) */
    double nu;         /* Vol of vol */
    double rho;        /* Correlation */
} fdp_sabr_params_t;

fdp_model_t* fdp_model_new_sabr_internal(
    fdp_context_t* ctx,
    double rate,
    double forward,
    double alpha,
    double beta,
    double nu,
    double rho
);

#endif /* FDPRICING_INTERNAL_MODEL_H */
