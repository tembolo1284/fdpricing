/* include/internal/utils/allocator.h - Memory allocation utilities */

#ifndef FDP_INTERNAL_ALLOCATOR_H
#define FDP_INTERNAL_ALLOCATOR_H

#include "fdpricing.h"

/* No need to redefine fdp_context_t - it's already in fdpricing.h */

/* ========================================================================
 * Context-aware allocation macros
 * 
 * These are defined as macros to avoid circular dependencies.
 * The actual inline functions are in context.h after the struct is defined.
 * ======================================================================== */

/* Note: Use fdp_ctx_malloc, fdp_ctx_realloc, fdp_ctx_free functions
 * defined in context.h for actual implementation */

/* Convenience macro for allocating arrays */
#define FDP_ALLOC_ARRAY(type, count) \
    ((type*)fdp_malloc(sizeof(type) * (count)))

/* Global allocation wrappers - these call the public API */
#define FDP_MALLOC(size)           fdp_malloc(size)
#define FDP_REALLOC(ptr, size)     fdp_realloc(ptr, size)
#define FDP_CALLOC(count, size)    fdp_calloc(count, size)
#define FDP_FREE(ptr)              fdp_free(ptr)

#endif /* FDP_INTERNAL_ALLOCATOR_H */
