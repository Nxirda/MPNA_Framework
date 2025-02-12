#pragma once 

#include "matrices_utils.h"
#include "vector.h"

typedef struct csr_matrix_s
{
    usz size;
    usz dim_x;
    usz dim_y;

    usz *col_index;
    usz *row_index; 
    
    f64 *data;
} csr_matrix_t;

//extern usz count_NNZ_elements(usz mesh_size, usz matrix_size);
//extern void allocate_CSR(usz mesh_size, csr_matrix_t *matrix);
//extern void fill_CSR(usz mesh_size, csr_matrix_t *matrix);
extern void allocate_CSR(usz dim_x, usz dim_y, usz NNZ, csr_matrix_t *matrix);
extern void print_CSR(csr_matrix_t *matrix);
extern void free_CSR(csr_matrix_t *matrix);
