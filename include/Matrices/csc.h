#pragma once

#include "matrices_utils.h"

typedef struct csc_matrix_s
{
    usz size;
    usz *row_index;
    usz *col_index;
    f64 *data;
} csc_matrix_t;

extern usz count_NNZ_elements_2(usz mesh_size, usz matrix_size);
extern void allocate_CSC(usz mesh_size, csc_matrix_t *matrix);
extern void fill_CSC(usz mesh_size, csc_matrix_t *matrix);
extern void print_CSC(csc_matrix_t *matrix);
extern void free_CSC(csc_matrix_t *matrix);
