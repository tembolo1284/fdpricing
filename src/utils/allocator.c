/* src/utils/allocator.c - Global memory allocator implementation */

#include "fdpricing.h"
#include <stdlib.h>
#include <string.h>

/* ========================================================================
 * Global Allocators (can be overridden)
 * ======================================================================== */

static struct {
    void* (*f_malloc)(size_t);
    void* (*f_realloc)(void*, size_t);
    void  (*f_free)(void*);
    int initialized;
} g_allocators = {
    .f_malloc = NULL,
    .f_realloc = NULL,
    .f_free = NULL,
    .initialized = 0
};

/* ========================================================================
 * Set Global Allocators
 * ======================================================================== */

void fdp_set_allocators(
    void* (*f_malloc)(size_t),
    void* (*f_realloc)(void*, size_t),
    void (*f_free)(void*))
{
    if (f_malloc && f_realloc && f_free) {
        g_allocators.f_malloc = f_malloc;
        g_allocators.f_realloc = f_realloc;
        g_allocators.f_free = f_free;
        g_allocators.initialized = 1;
    } else {
        /* Reset to default */
        g_allocators.f_malloc = NULL;
        g_allocators.f_realloc = NULL;
        g_allocators.f_free = NULL;
        g_allocators.initialized = 0;
    }
}

/* ========================================================================
 * Global Memory Functions
 * ======================================================================== */

void* fdp_malloc(size_t size)
{
    if (g_allocators.initialized && g_allocators.f_malloc) {
        return g_allocators.f_malloc(size);
    }
    return malloc(size);
}

void* fdp_realloc(void* ptr, size_t size)
{
    if (g_allocators.initialized && g_allocators.f_realloc) {
        return g_allocators.f_realloc(ptr, size);
    }
    return realloc(ptr, size);
}

void* fdp_calloc(size_t count, size_t size)
{
    if (g_allocators.initialized && g_allocators.f_malloc) {
        size_t total = count * size;
        void* ptr = g_allocators.f_malloc(total);
        if (ptr) {
            memset(ptr, 0, total);
        }
        return ptr;
    }
    return calloc(count, size);
}

void fdp_free(void* ptr)
{
    if (!ptr) return;
    
    if (g_allocators.initialized && g_allocators.f_free) {
        g_allocators.f_free(ptr);
    } else {
        free(ptr);
    }
}
