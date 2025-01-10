#pragma once 

#include "matrices_utils.h"

typedef struct csr_matrix_s
{
    usz size;
    usz *col_index;
    usz *row_index; 
    f64 *values;
} csr_matrix_t;

extern usz count_NNZ_elements(usz mesh_size, usz matrix_size);
extern void allocate_CSR(usz mesh_size, csr_matrix_t *matrix);
extern void fill_CSR(usz mesh_size, csr_matrix_t *matrix);
extern void print_CSR(csr_matrix_t *matrix);
extern void free_CSR(csr_matrix_t *matrix);
