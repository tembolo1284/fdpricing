/**
 * context.c - Context management implementation
 */

#include "fdpricing.h"
#include "internal/core/context.h"
#include "internal/utils/allocator.h"

fdp_context_t* fdp_context_new(void)
{
    fdp_context_t* ctx = FDP_ALLOC_ARRAY(fdp_context_t, 1);
    if (!ctx) {
        return NULL;
    }
    
    /* Initialize with defaults */
    ctx->allocators.f_malloc = NULL;
    ctx->allocators.f_realloc = NULL;
    ctx->allocators.f_free = NULL;
    ctx->use_custom_allocators = 0;
    ctx->last_error = FDP_SUCCESS;
    
    return ctx;
}

void fdp_context_free(fdp_context_t* ctx)
{
    if (!ctx) return;
    
    /* Note: We always free the context itself using global allocator,
     * since it was allocated with global allocator */
    fdp_free(ctx);
}

void fdp_context_set_allocators(
    fdp_context_t* ctx,
    void* (*f_malloc)(size_t),
    void* (*f_realloc)(void*, size_t),
    void (*f_free)(void*))
{
    if (!ctx) return;
    
    if (f_malloc && f_realloc && f_free) {
        ctx->allocators.f_malloc = f_malloc;
        ctx->allocators.f_realloc = f_realloc;
        ctx->allocators.f_free = f_free;
        ctx->use_custom_allocators = 1;
    } else {
        /* Reset to global allocators */
        ctx->allocators.f_malloc = NULL;
        ctx->allocators.f_realloc = NULL;
        ctx->allocators.f_free = NULL;
        ctx->use_custom_allocators = 0;
    }
}

/* Error handling */
fdp_error_t fdp_get_last_error(fdp_context_t* ctx)
{
    return fdp_ctx_get_error(ctx);
}

void fdp_clear_error(fdp_context_t* ctx)
{
    if (ctx) {
        ctx->last_error = FDP_SUCCESS;
    }
}

const char* fdp_get_error_string(fdp_error_t error)
{
    switch (error) {
        case FDP_SUCCESS:
            return "Success";
        case FDP_ERROR_INVALID_PARAM:
            return "Invalid parameter";
        case FDP_ERROR_ALLOCATION:
            return "Memory allocation failed";
        case FDP_ERROR_CONVERGENCE:
            return "Solver failed to converge";
        case FDP_ERROR_NOT_IMPLEMENTED:
            return "Feature not yet implemented";
        case FDP_ERROR_STABILITY:
            return "Numerical stability violation (check CFL condition)";
        default:
            return "Unknown error";
    }
}
