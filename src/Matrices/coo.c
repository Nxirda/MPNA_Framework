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
    printf("Not done\n");
    return;
}
/*assert(matrix->data != NULL);
 
    usz idx = 0;

    for(usz i = 0; i < matrix->dim_x; i++)
    {
        usz elems_in_row = matrix->row_index[i+1] - matrix->row_index[i]; 
        
        for(usz j = 0; j < matrix->dim_y; j++)
        {
            if(elems_in_row > 0 && matrix->col_index[idx] == j)
            {
                printf("%-2.3f ", matrix->data[idx]);
                elems_in_row --;
                idx ++;
            }
            else
            {
                printf("__.___ " );
        
            }
        }
        printf("\n");

    }
}*/

void free_COO(coo_matrix_t *matrix)
{
    assert(matrix != NULL);

    free(matrix->col_index);
    free(matrix->row_index);
    free(matrix->data);
}
