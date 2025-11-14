/**
 * fdpricing.h - Finite Difference Option Pricing Library
 * 
 * A pure C library for pricing options using finite difference methods.
 * Supports multiple models (GBM, Heston, SABR, Jump-Diffusion) and
 * option types (European, American, Bermudan, Barrier, Asian, Lookback).
 */

#ifndef FDPRICING_H_INCLUDED
#define FDPRICING_H_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* Version information */
#define FDP_VERSION_MAJOR 0
#define FDP_VERSION_MINOR 1
#define FDP_VERSION_PATCH 0
#define FDP_VERSION ((FDP_VERSION_MAJOR << 16) | (FDP_VERSION_MINOR << 8) | FDP_VERSION_PATCH)

/* Export macros for shared library */
#ifndef FDP_API
# ifdef _WIN32
#  if defined(FDP_BUILD_SHARED)
#   define FDP_API __declspec(dllexport)
#  elif !defined(FDP_BUILD_STATIC)
#   define FDP_API __declspec(dllimport)
#  else
#   define FDP_API
#  endif
# else
#  if __GNUC__ >= 4
#   define FDP_API __attribute__((visibility("default")))
#  else
#   define FDP_API
#  endif
# endif
#endif

/* ========================================================================
 * Version and Compatibility
 * ======================================================================== */

FDP_API uint32_t fdp_get_version(void);
FDP_API int fdp_is_compatible(void);

/* ========================================================================
 * Memory Allocators
 * ======================================================================== */

FDP_API void fdp_set_allocators(
    void* (*f_malloc)(size_t),
    void* (*f_realloc)(void*, size_t),
    void (*f_free)(void*)
);

FDP_API void* fdp_malloc(size_t size);
FDP_API void* fdp_realloc(void* ptr, size_t size);
FDP_API void* fdp_calloc(size_t count, size_t size);
FDP_API void fdp_free(void* ptr);

/* ========================================================================
 * Core Types (Opaque)
 * ======================================================================== */

struct fdp_context_s;
typedef struct fdp_context_s fdp_context_t;

struct fdp_grid_s;
typedef struct fdp_grid_s fdp_grid_t;

struct fdp_model_s;
typedef struct fdp_model_s fdp_model_t;

struct fdp_option_s;
typedef struct fdp_option_s fdp_option_t;

struct fdp_solver_params_s;
typedef struct fdp_solver_params_s fdp_solver_params_t;

struct fdp_result_s;
typedef struct fdp_result_s fdp_result_t;

/* ========================================================================
 * Enumerations
 * ======================================================================== */

/* Option types */
typedef enum {
    FDP_OPTION_CALL = 0,
    FDP_OPTION_PUT  = 1
} fdp_option_type_t;

/* Option styles */
typedef enum {
    FDP_STYLE_EUROPEAN = 0,
    FDP_STYLE_AMERICAN = 1,
    FDP_STYLE_BERMUDAN = 2
} fdp_option_style_t;

/* Barrier types */
typedef enum {
    FDP_BARRIER_NONE      = 0,
    FDP_BARRIER_UP_OUT    = 1,
    FDP_BARRIER_UP_IN     = 2,
    FDP_BARRIER_DOWN_OUT  = 3,
    FDP_BARRIER_DOWN_IN   = 4
} fdp_barrier_type_t;

/* Model types */
typedef enum {
    FDP_MODEL_GBM          = 0,
    FDP_MODEL_HESTON       = 1,
    FDP_MODEL_SABR         = 2,
    FDP_MODEL_MERTON_JUMP  = 3,
    FDP_MODEL_KOU_JUMP     = 4,
    FDP_MODEL_BATES        = 5,
    FDP_MODEL_LOCAL_VOL    = 6
} fdp_model_type_t;

/* Solver methods */
typedef enum {
    FDP_SOLVER_EXPLICIT       = 0, /* Explicit Euler (forward)        */
    FDP_SOLVER_IMPLICIT       = 1, /* Implicit Euler (backward)       */
    FDP_SOLVER_CRANK_NICOLSON = 2, /* Crank-Nicolson                  */
    FDP_SOLVER_PSOR           = 3  /* Projected SOR (for American)    */
} fdp_solver_method_t;

/* Grid types */
typedef enum {
    FDP_GRID_UNIFORM      = 0,
    FDP_GRID_LOG_UNIFORM  = 1,
    FDP_GRID_SINH         = 2,
    FDP_GRID_ADAPTIVE     = 3
} fdp_grid_type_t;

/* Error codes */
typedef enum {
    FDP_SUCCESS           = 0,
    FDP_ERROR_INVALID_PARAM   = -1,
    FDP_ERROR_ALLOCATION      = -2,
    FDP_ERROR_CONVERGENCE     = -3,
    FDP_ERROR_NOT_IMPLEMENTED = -4,
    FDP_ERROR_STABILITY       = -5
} fdp_error_t;

/* ========================================================================
 * Context Management
 * ======================================================================== */

FDP_API fdp_context_t* fdp_context_new(void);
FDP_API void fdp_context_free(fdp_context_t* ctx);
FDP_API void fdp_context_set_allocators(
    fdp_context_t* ctx,
    void* (*f_malloc)(size_t),
    void* (*f_realloc)(void*, size_t),
    void (*f_free)(void*)
);

/* ========================================================================
 * Grid Management
 * ======================================================================== */

/**
 * Create a 1D grid for spatial discretization
 * 
 * @param ctx Context object
 * @param grid_type Type of grid spacing
 * @param n_space Number of spatial points
 * @param n_time Number of time steps
 * @param spot Current spot price
 * @param s_min Minimum S value
 * @param s_max Maximum S value
 * @param t_max Time to maturity
 */
FDP_API fdp_grid_t* fdp_grid_new_1d(
    fdp_context_t* ctx,
    fdp_grid_type_t grid_type,
    int n_space,
    int n_time,
    double spot,
    double s_min,
    double s_max,
    double t_max
);

/**
 * Create a 2D grid (for stochastic volatility models)
 */
FDP_API fdp_grid_t* fdp_grid_new_2d(
    fdp_context_t* ctx,
    int n_space,
    int n_vol,
    int n_time,
    double spot,
    double s_min,
    double s_max,
    double vol,
    double v_min,
    double v_max,
    double t_max
);

FDP_API void fdp_grid_free(fdp_grid_t* grid);

/* Grid accessors */
FDP_API int fdp_grid_get_n_space(const fdp_grid_t* grid);
FDP_API int fdp_grid_get_n_time(const fdp_grid_t* grid);
FDP_API double fdp_grid_get_spot(const fdp_grid_t* grid);
FDP_API const double* fdp_grid_get_space_points(const fdp_grid_t* grid);
FDP_API const double* fdp_grid_get_time_points(const fdp_grid_t* grid);

/* ========================================================================
 * Model Management
 * ======================================================================== */

/* GBM (Black-Scholes) Model */
FDP_API fdp_model_t* fdp_model_new_gbm(
    fdp_context_t* ctx,
    double rate,      /* Risk-free rate   */
    double div_yield, /* Dividend yield   */
    double vol        /* Volatility       */
);

/* Heston Stochastic Volatility Model */
FDP_API fdp_model_t* fdp_model_new_heston(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double kappa,     /* Mean reversion speed */
    double theta,     /* Long-term variance   */
    double sigma,     /* Vol of vol           */
    double rho,       /* Correlation          */
    double v0         /* Initial variance     */
);

/* Merton Jump-Diffusion Model */
FDP_API fdp_model_t* fdp_model_new_merton_jump(
    fdp_context_t* ctx,
    double rate,
    double div_yield,
    double vol,
    double lambda,    /* Jump intensity      */
    double mu_j,      /* Jump mean (log)     */
    double sigma_j    /* Jump std dev (log)  */
);

FDP_API void fdp_model_free(fdp_model_t* model);
FDP_API fdp_model_type_t fdp_model_get_type(const fdp_model_t* model);

/* ========================================================================
 * Option Specification
 * ======================================================================== */

/* Vanilla European/American option */
FDP_API fdp_option_t* fdp_option_new_vanilla(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    fdp_option_style_t style,
    double strike,
    double maturity
);

/* Bermudan option (discrete exercise dates) */
FDP_API fdp_option_t* fdp_option_new_bermudan(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,
    double maturity,
    const double* exercise_times,
    int n_exercise_times
);

/* Barrier option */
FDP_API fdp_option_t* fdp_option_new_barrier(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,
    double maturity,
    fdp_barrier_type_t barrier_type,
    double barrier_level,
    double rebate
);

/* Asian option (arithmetic average) */
FDP_API fdp_option_t* fdp_option_new_asian(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,
    double maturity,
    int n_averaging_points
);

/* Lookback option */
FDP_API fdp_option_t* fdp_option_new_lookback(
    fdp_context_t* ctx,
    fdp_option_type_t type,
    double strike,     /* For fixed strike, 0.0 for floating */
    double maturity,
    int is_fixed_strike
);

FDP_API void fdp_option_free(fdp_option_t* option);

/* ========================================================================
 * Solver Parameters
 * ======================================================================== */

FDP_API fdp_solver_params_t* fdp_solver_params_new(fdp_context_t* ctx);
FDP_API void fdp_solver_params_free(fdp_solver_params_t* params);

FDP_API void fdp_solver_params_set_method(
    fdp_solver_params_t* params,
    fdp_solver_method_t method
);

FDP_API void fdp_solver_params_set_theta(
    fdp_solver_params_t* params,
    double theta  /* 0.0 = explicit, 0.5 = CN, 1.0 = implicit */
);

FDP_API void fdp_solver_params_set_tolerance(
    fdp_solver_params_t* params,
    double tolerance
);

FDP_API void fdp_solver_params_set_max_iterations(
    fdp_solver_params_t* params,
    int max_iter
);

FDP_API void fdp_solver_params_set_omega(
    fdp_solver_params_t* params,
    double omega  /* SOR relaxation parameter */
);

/* ========================================================================
 * PDE Solver (Main Interface)
 * ======================================================================== */

/**
 * Solve the PDE for the given option, model, and grid
 * 
 * Returns a result object containing the solution surface
 */
FDP_API fdp_result_t* fdp_solve_pde(
    fdp_context_t* ctx,
    const fdp_model_t* model,
    const fdp_option_t* option,
    const fdp_grid_t* grid,
    const fdp_solver_params_t* params
);

FDP_API void fdp_result_free(fdp_result_t* result);

/* Result accessors */
FDP_API double fdp_result_get_price(const fdp_result_t* result, double spot);
FDP_API double fdp_result_get_delta(const fdp_result_t* result, double spot);
FDP_API double fdp_result_get_gamma(const fdp_result_t* result, double spot);
FDP_API double fdp_result_get_theta(const fdp_result_t* result, double spot);
FDP_API int    fdp_result_get_iterations(const fdp_result_t* result);
FDP_API fdp_error_t fdp_result_get_error_code(const fdp_result_t* result);

/* Get the full solution surface (n_space × n_time array) */
FDP_API const double* fdp_result_get_surface(const fdp_result_t* result);

/* ========================================================================
 * Convenience Functions (Simple Interface)
 * ======================================================================== */

/**
 * Simple interface for common use cases
 * Returns price directly, no need to manage objects
 */

/* European vanilla */
FDP_API double fdp_price_european_call(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

FDP_API double fdp_price_european_put(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

/* American vanilla */
FDP_API double fdp_price_american_call(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

FDP_API double fdp_price_american_put(
    double spot,
    double strike,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

/* Barrier options: generic convenience function */
FDP_API double fdp_price_barrier_option(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    fdp_option_type_t option_type,
    fdp_barrier_type_t barrier_type,
    double rebate,
    int n_space,
    int n_time
);

/* Barrier options: specific wrappers */

/* Up-and-out call (no rebate) */
FDP_API double fdp_price_up_and_out_call(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

/* Down-and-out put (no rebate) */
FDP_API double fdp_price_down_and_out_put(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

/* Up-and-in call (no rebate) */
FDP_API double fdp_price_up_and_in_call(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

/* Down-and-in put (no rebate) */
FDP_API double fdp_price_down_and_in_put(
    double spot,
    double strike,
    double barrier,
    double rate,
    double div_yield,
    double vol,
    double maturity,
    int n_space,
    int n_time
);

/* ========================================================================
 * Error Handling
 * ======================================================================== */

FDP_API const char*   fdp_get_error_string(fdp_error_t error);
FDP_API fdp_error_t   fdp_get_last_error(fdp_context_t* ctx);
FDP_API void          fdp_clear_error(fdp_context_t* ctx);

#endif /* FDPRICING_H_INCLUDED */

