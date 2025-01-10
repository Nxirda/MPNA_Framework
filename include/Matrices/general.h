#pragma once

#include "matrices_utils.h"

typedef struct matrix_s
{
    usz size;
    f64 *data;
} matrix_t;

extern void allocate_matrix(usz mesh_size, matrix_t *matrix);
extern void fill_matrix(usz mesh_size, matrix_t matrix);
extern void print_matrix(matrix_t matrix);
extern void free_matrix(matrix_t matrix);
