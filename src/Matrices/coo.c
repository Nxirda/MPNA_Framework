#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "coo.h"

void allocate_COO(usz dim_x, usz dim_y, usz NNZ, coo_matrix_t *matrix)
{
    const usz N = dim_x * dim_y;
    matrix->size = N;

    matrix->dim_x = dim_x;
    matrix->dim_y = dim_y;

    matrix->data        = (f64 *)malloc(NNZ * sizeof(f64));
    matrix->col_index   = (usz *)malloc(NNZ * sizeof(usz));
    matrix->row_index   = (usz *)malloc(NNZ * sizeof(usz));
}


void print_COO(coo_matrix_t *matrix)
{
    assert(matrix->data != NULL);
 
    usz size = sizeof(matrix->data)/sizeof(matrix->data[0]);
    printf("SIZE = %d\n", size);
    for(usz i = 0; i < size; i++)
    {
        printf("row: %zu col: %zu value: %f\n", matrix->row_index[i], 
                matrix->col_index[i], matrix->data[i]);
    }
}

void free_COO(coo_matrix_t *matrix)
{
    assert(matrix != NULL);

    free(matrix->col_index);
    free(matrix->row_index);
    free(matrix->data);
}
