/* include/internal/core/context.h - Context structure and utilities */

#ifndef FDP_INTERNAL_CONTEXT_H
#define FDP_INTERNAL_CONTEXT_H

#include "fdpricing.h"

/* ========================================================================
 * Function Pointer Types
 * ======================================================================== */

typedef void* (*fdp_malloc_fn)(size_t);
typedef void* (*fdp_realloc_fn)(void*, size_t);
typedef void  (*fdp_free_fn)(void*);

/* ========================================================================
 * Context Structure Definition
 * ======================================================================== */

struct fdp_context_s {
    /* Custom allocators */
    struct {
        fdp_malloc_fn f_malloc;
        fdp_realloc_fn f_realloc;
        fdp_free_fn f_free;
    } allocators;
    
    int use_custom_allocators;
    
    /* Error state */
    fdp_error_t last_error;
};

/* ========================================================================
 * Context-aware Allocation Functions
 * ======================================================================== */

static inline void* fdp_ctx_malloc(fdp_context_t* ctx, size_t size)
{
    if (ctx && ctx->use_custom_allocators && ctx->allocators.f_malloc) {
        return ctx->allocators.f_malloc(size);
    }
    return fdp_malloc(size);
}

static inline void* fdp_ctx_realloc(fdp_context_t* ctx, void* ptr, size_t size)
{
    if (ctx && ctx->use_custom_allocators && ctx->allocators.f_realloc) {
        return ctx->allocators.f_realloc(ptr, size);
    }
    return fdp_realloc(ptr, size);
}

static inline void fdp_ctx_free(fdp_context_t* ctx, void* ptr)
{
    if (!ptr) return;
    
    if (ctx && ctx->use_custom_allocators && ctx->allocators.f_free) {
        ctx->allocators.f_free(ptr);
    } else {
        fdp_free(ptr);
    }
}

/* ========================================================================
 * Context-aware Allocation Macros
 * ======================================================================== */

#define FDP_CTX_ALLOC(ctx, size) \
    fdp_ctx_malloc(ctx, size)

#define FDP_CTX_REALLOC(ctx, ptr, size) \
    fdp_ctx_realloc(ctx, ptr, size)

#define FDP_CTX_FREE(ctx, ptr) \
    fdp_ctx_free(ctx, ptr)

#define FDP_CTX_ALLOC_ARRAY(ctx, type, count) \
    ((type*)fdp_ctx_malloc(ctx, sizeof(type) * (count)))

#define FDP_CTX_CALLOC_ARRAY(ctx, type, count) \
    ((type*)fdp_calloc(count, sizeof(type)))

/* ========================================================================
 * Error Management
 * ======================================================================== */

static inline void fdp_ctx_set_error(fdp_context_t* ctx, fdp_error_t error)
{
    if (ctx) {
        ctx->last_error = error;
    }
}

static inline fdp_error_t fdp_ctx_get_error(const fdp_context_t* ctx)
{
    return ctx ? ctx->last_error : FDP_ERROR_INVALID_PARAM;
}

#endif /* FDP_INTERNAL_CONTEXT_H */
