#pragma once

#include "matrices_utils.h"

typedef struct matrix_s
{
    usz size;
    usz dim_x;
    usz dim_y;

    f64 *data;
} matrix_t;

extern void allocate_matrix(usz size_x, usz size_y, matrix_t *matrix);
//extern void allocate_matrix(usz mesh_size, matrix_t *matrix);
//extern void fill_matrix(usz mesh_size, matrix_t matrix);
extern void print_matrix(matrix_t matrix);
extern void free_matrix(matrix_t matrix);
