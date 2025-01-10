#pragma once

#include "types.h"

typedef struct vector_s
{
    usz size;
    f64 *data;
} vector_t;

extern void allocate_vector(vector_t *vector, usz size);
extern void fill_vector(vector_t *vector, usz begin);
extern void print_vector(vector_t *vector);
extern void free_vector(vector_t *vector);
extern void swap_vector(vector_t *a, vector_t *b);
