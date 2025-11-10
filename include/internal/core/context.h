/**
 * context.h - Context object internal structure
 * 
 * The context holds all per-instance state for the library.
 * This allows multiple independent contexts without global state.
 */

#ifndef FDPRICING_INTERNAL_CONTEXT_H
#define FDPRICING_INTERNAL_CONTEXT_H

#include "fdpricing.h"
#include "internal/utils/allocator.h"

/* Internal context structure */
struct fdp_context_s {
    /* Per-context allocators (optional, can be NULL to use global) */
    struct {
        fdp_malloc_fn f_malloc;
        fdp_realloc_fn f_realloc;
        fdp_free_fn f_free;
    } allocators;
    
    /* Error state */
    fdp_error_t last_error;
    
    /* Configuration flags */
    int use_custom_allocators;  /* 0 = use global, 1 = use context-specific */
};

/* Context-aware allocation functions */
static inline void* fdp_ctx_malloc(fdp_context_t* ctx, size_t size)
{
    if (ctx && ctx->use_custom_allocators) {
        return ctx->allocators.f_malloc(size);
    }
    return fdp_malloc(size);
}

static inline void* fdp_ctx_realloc(fdp_context_t* ctx, void* ptr, size_t size)
{
    if (ctx && ctx->use_custom_allocators) {
        return ctx->allocators.f_realloc(ptr, size);
    }
    return fdp_realloc(ptr, size);
}

static inline void* fdp_ctx_calloc(fdp_context_t* ctx, size_t count, size_t size)
{
    size_t total = count * size;
    void* ptr = fdp_ctx_malloc(ctx, total);
    if (ptr) {
        char* p = (char*)ptr;
        for (size_t i = 0; i < total; ++i) {
            p[i] = 0;
        }
    }
    return ptr;
}

static inline void fdp_ctx_free(fdp_context_t* ctx, void* ptr)
{
    if (!ptr) return;
    
    if (ctx && ctx->use_custom_allocators) {
        ctx->allocators.f_free(ptr);
    } else {
        fdp_free(ptr);
    }
}

/* Array allocation with context */
#define FDP_CTX_ALLOC_ARRAY(ctx, type, count) \
    ((type*)fdp_ctx_malloc((ctx), (count) * sizeof(type)))

#define FDP_CTX_CALLOC_ARRAY(ctx, type, count) \
    ((type*)fdp_ctx_calloc((ctx), (count), sizeof(type)))

/* Error handling helpers */
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

#endif /* FDPRICING_INTERNAL_CONTEXT_H */
