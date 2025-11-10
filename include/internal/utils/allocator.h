/**
 * allocator.h - Internal memory allocation interface
 * 
 * Provides a global allocator system that can be customized by the user.
 * All internal code should use these functions instead of malloc/free directly.
 */

#ifndef FDPRICING_INTERNAL_ALLOCATOR_H
#define FDPRICING_INTERNAL_ALLOCATOR_H

#include <stddef.h>

/* Allocator function pointer types */
typedef void* (*fdp_malloc_fn)(size_t);
typedef void* (*fdp_realloc_fn)(void*, size_t);
typedef void (*fdp_free_fn)(void*);

/* Global allocator state */
struct fdp_allocators_s {
    fdp_malloc_fn f_malloc;
    fdp_realloc_fn f_realloc;
    fdp_free_fn f_free;
};

/* Declared here, defined in allocator.c */
extern struct fdp_allocators_s _fdp_allocators;

/* Internal allocation functions (inline for performance) */
static inline void* fdp_malloc(size_t size)
{
    return _fdp_allocators.f_malloc(size);
}

static inline void* fdp_realloc(void* ptr, size_t size)
{
    return _fdp_allocators.f_realloc(ptr, size);
}

static inline void fdp_free(void* ptr)
{
    _fdp_allocators.f_free(ptr);
}

/* Utility allocation functions */
static inline void* fdp_calloc(size_t count, size_t size)
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

static inline char* fdp_strdup(const char* str)
{
    if (!str) return NULL;
    
    size_t len = 0;
    while (str[len]) ++len;
    len++; /* Include null terminator */
    
    char* result = (char*)fdp_malloc(len);
    if (result) {
        for (size_t i = 0; i < len; ++i) {
            result[i] = str[i];
        }
    }
    return result;
}

/* Array allocation helpers with type safety */
#define FDP_ALLOC_ARRAY(type, count) \
    ((type*)fdp_malloc((count) * sizeof(type)))

#define FDP_CALLOC_ARRAY(type, count) \
    ((type*)fdp_calloc((count), sizeof(type)))

#define FDP_REALLOC_ARRAY(ptr, type, count) \
    ((type*)fdp_realloc((ptr), (count) * sizeof(type)))

/* Initialize allocators (called internally) */
void fdp_allocators_init(void);

/* Set custom allocators (implementation of public API) */
void fdp_allocators_set(
    fdp_malloc_fn f_malloc,
    fdp_realloc_fn f_realloc,
    fdp_free_fn f_free
);

#endif /* FDPRICING_INTERNAL_ALLOCATOR_H */
