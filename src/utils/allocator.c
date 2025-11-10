/**
 * allocator.c - Memory allocator implementation
 */

#include "fdpricing.h"
#include "internal/utils/allocator.h"
#include <stdlib.h>

/* Global allocator state - initialized with standard library functions */
struct fdp_allocators_s _fdp_allocators = {
    malloc,
    realloc,
    free
};

void fdp_allocators_init(void)
{
    /* Reset to default allocators */
    _fdp_allocators.f_malloc = malloc;
    _fdp_allocators.f_realloc = realloc;
    _fdp_allocators.f_free = free;
}

void fdp_allocators_set(
    fdp_malloc_fn f_malloc,
    fdp_realloc_fn f_realloc,
    fdp_free_fn f_free)
{
    if (f_malloc && f_realloc && f_free) {
        _fdp_allocators.f_malloc = f_malloc;
        _fdp_allocators.f_realloc = f_realloc;
        _fdp_allocators.f_free = f_free;
    }
}

/* ========================================================================
 * Public API Implementation
 * ======================================================================== */

void fdp_set_allocators(
    void* (*f_malloc)(size_t),
    void* (*f_realloc)(void*, size_t),
    void (*f_free)(void*))
{
    fdp_allocators_set(f_malloc, f_realloc, f_free);
}

/* These just forward to the inline functions, but need to exist for linking */
void* fdp_malloc(size_t size)
{
    return _fdp_allocators.f_malloc(size);
}

void* fdp_realloc(void* ptr, size_t size)
{
    return _fdp_allocators.f_realloc(ptr, size);
}

void* fdp_calloc(size_t count, size_t size)
{
    size_t total = count * size;
    void* ptr = fdp_malloc(total);
    if (ptr) {
        char* p = (char*)ptr;
        for (size_t i = 0; i < total; ++i) {
            p[i] = 0;
        }
    }
    return ptr;
}

void fdp_free(void* ptr)
{
    if (ptr) {
        _fdp_allocators.f_free(ptr);
    }
}
