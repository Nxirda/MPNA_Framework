#pragma once

#include "types.h"

typedef struct vector_s
{
    usz size;
    f64 *data;
} vector_t;

static inline f64 f64_random(f64 min, f64 max);
static inline f64 f64_abs(f64 num);

// Allocation function
extern void allocate_vector(vector_t *vector, usz size);
// Might add ine that enforces alignment on specific boundary

// Intializations functions
extern void init_iota_vector(vector_t *vector, f64 begin);
extern void init_random_vector(vector_t *a, f64 min, f64 max);
extern void init_constant_vector(vector_t *a, f64 value);

// Utilities functions
extern void copy_vector(vector_t const*src, vector_t *dest);
extern void print_vector(vector_t *vector);
extern void free_vector(vector_t *vector);
extern void swap_vector(vector_t *a, vector_t *b);

// Maths functions
extern f64 dot_product(vector_t *a, vector_t *b);
extern void add_vector(vector_t *a, vector_t *b, vector_t *c);
extern void mul_vector(vector_t *a, vector_t *b, vector_t *c);

// Scalar x vector operations
extern void add_scalar_vector(f64 scalar, vector_t *b, vector_t *c);
extern void mul_scalar_vector(f64 scalar, vector_t *b, vector_t *c);
extern void daxpy(f64 scalar, vector_t *a, vector_t *b, vector_t *c);

extern u8 equal_vector(vector_t *a, vector_t *b);
extern u8 allclose_vector(vector_t *a, vector_t *b, f64 atol, f64 rtol);
