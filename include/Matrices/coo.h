#pragma once

#include "matrices_utils.h"

typedef struct coo_matrix_s
{
    usz size;
    usz dim_x;
    usz dim_y;

    usz *col_index;
    usz *row_index; 
    
    f64 *data;
} coo_matrix_t;

extern void allocate_COO(usz dim_x, usz dim_y, usz NNZ, coo_matrix_t *matrix);
extern void print_COO(coo_matrix_t *matrix);
extern void free_COO(coo_matrix_t *matrix);
