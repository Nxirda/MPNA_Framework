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
extern void init_random_vector(vector_t *vector, f64 min, f64 max);
extern void init_constant_vector(vector_t *vector, f64 value);

// Utilities functions
extern void copy_vector(vector_t const *src, vector_t *dest);
extern void print_vector(vector_t const *vector);
extern void free_vector(vector_t *vector);
extern void swap_vector(vector_t *a, vector_t *b);

// Maths functions (Equivalent to BLAS 1 operations)
extern f64 dot_product(vector_t const *a, vector_t const *b);
extern f64 dot_product_conjugate(vector_t const *a, vector_t const *b);
extern f64 norm_one_vector();
extern f64 norm_two_vector();

extern void rotate_vector();
// See rot, rotg, rotm, rotmg 
extern void add_vector(vector_t const *a, vector_t const *b, vector_t *dest);
extern void mul_vector(vector_t const *a, vector_t const *b, vector_t *dest);

extern void add_scalar_vector(f64 scalar, vector_t const *b, vector_t *dest);
extern void mul_scalar_vector(f64 scalar, vector_t const *b, vector_t *dest);
extern void daxpy(f64 scalar, vector_t const *a, vector_t const *b, vector_t *dest);

extern u8 equal_vector(vector_t const *a, vector_t const *b);
extern u8 allclose_vector(vector_t const *a, vector_t const *b, f64 atol, f64 rtol);
