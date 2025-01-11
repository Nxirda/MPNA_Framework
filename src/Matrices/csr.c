#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "csr.h"

usz count_NNZ_elements(usz mesh_size, usz matrix_size)
{
    assert(mesh_size > 0 && matrix_size > 0);

    const usz other_diagonal_size = matrix_size - mesh_size; 
    return matrix_size + 4 * other_diagonal_size; 
}

void allocate_CSR(usz mesh_size, csr_matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);

    const usz N = mesh_size * mesh_size;
    const usz NNZ = count_NNZ_elements(mesh_size, N);

    matrix->size = N;
    matrix->data        = (f64 *)malloc(NNZ * sizeof(f64));
    matrix->col_index   = (usz *)malloc(NNZ    * sizeof(usz));
    matrix->row_index   = (usz *)malloc((N+1)  * sizeof(usz));

    //usz bytes_data = NNZ * sizeof(f64);
    //usz bytes_col = NNZ * sizeof(usz);
    //usz bytes_row = (N+1) * sizeof(usz);
    //usz total = bytes_data + bytes_col + bytes_row;
    //printf("Allocated size : %ld, %ld, %ld, total : %ld\n", 
    //        bytes_data, bytes_col, bytes_row, total);
}

/*
 * Starting from a mesh of size N by N 
 * Returns the CSR representation as a N^2 * N^2 matrix
 * */
void fill_CSR(usz mesh_size, csr_matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);

    usz CSR_idx = 0;
    usz row_idx = 0;
    
    matrix->row_index[row_idx] = 0;
    row_idx ++;

    for(usz j = 0; j < mesh_size; j++)
    {
        for(usz i = 0; i < mesh_size; i++) 
        {
           //usz old_Idx = CSR_idx;

            //Far left 
            if(j > 0)
            {
                usz col = to_matrix_id(i,j-1, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }
            //Close left
            if(i > 0)
            {
                usz col = to_matrix_id(i-1,j, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }
           
            usz col = to_matrix_id(i, j, mesh_size);
            matrix->data[CSR_idx] = 4.0;
            matrix->col_index[CSR_idx] = col;
            CSR_idx ++;

            //Close right
            if(i < mesh_size -1)
            {
                usz col = to_matrix_id(i+1,j, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }
            //Far right
            if(j < mesh_size -1)
            {
                usz col = to_matrix_id(i,j+1, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }

            matrix->row_index[row_idx] = CSR_idx;
            row_idx ++;
        }
    }
}

void print_CSR(csr_matrix_t *matrix)
{
    assert(matrix->data != NULL);
 
    usz idx = 0;

    for(usz i = 0; i < matrix->size; i++)
    {
        usz elems_in_row = matrix->row_index[i+1] - matrix->row_index[i]; 
        
        for(usz j = 0; j < matrix->size; j++)
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
}

void free_CSR(csr_matrix_t *matrix)
{
    assert(matrix != NULL);

    free(matrix->col_index);
    free(matrix->row_index);
    free(matrix->data);
}
