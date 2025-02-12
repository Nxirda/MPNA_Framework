#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "csc.h"

/*usz count_NNZ_elements_2(usz mesh_size, usz matrix_size)
{
    assert(mesh_size > 0 && matrix_size > 0);

    const usz other_diagonal_size = matrix_size - mesh_size; 
    return matrix_size + 4 * other_diagonal_size; 
}
*/
/*void allocate_CSC(usz mesh_size, csc_matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);

    const usz N = mesh_size * mesh_size;
    const usz NNZ = count_NNZ_elements_2(mesh_size, N);
    
    matrix->size = N;
    matrix->data        = (f64 *)malloc(NNZ * sizeof(f64));
    matrix->row_index   = (usz *)malloc(NNZ    * sizeof(usz));
    matrix->col_index   = (usz *)malloc((N+1)  * sizeof(usz));

    //usz bytes_data = NNZ * sizeof(f64);
    //usz bytes_col = NNZ * sizeof(usz);
    //usz bytes_row = (N+1) * sizeof(usz);
    //usz total = bytes_data + bytes_col + bytes_row;
    //printf("Allocated size : %ld, %ld, %ld, total : %ld\n", 
    //        bytes_data, bytes_col, bytes_row, total);
}*/

void allocate_CSC(usz dim_x, usz dim_y, usz NNZ, csc_matrix_t *matrix)
{
    const usz N = dim_x * dim_y;
    matrix->size = N;

    matrix->data = (f64 *)malloc(NNZ * sizeof(f64));
    matrix->col_index = (usz *)malloc((dim_y+1) * sizeof(usz));
    matrix->row_index = (usz *)malloc(NNZ * sizeof(usz));
}

/*
 * Starting from a mesh of size N by N 
 * Returns the CSC representation as a N^2 * N^2 matrix
 * */
/*
void fill_CSC(usz mesh_size, csc_matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);

    usz CSC_idx = 0;
    usz col_idx = 0;
    
    matrix->col_index[col_idx] = 0;
    col_idx ++;

    for(usz j = 0; j < mesh_size; j++)
    {
        for(usz i = 0; i < mesh_size; i++) 
        {
            //Far left 
            if(j > 0)
            {
                usz row = to_matrix_id(i,j-1, mesh_size);
                matrix->data[CSC_idx] = -1.0;
                matrix->row_index[CSC_idx] = row;
                CSC_idx ++;
            }
            //Close left
            if(i > 0)
            {
                usz row = to_matrix_id(i-1,j, mesh_size);
                matrix->data[CSC_idx] = -1.0;
                matrix->row_index[CSC_idx] = row;
                CSC_idx ++;
            }
           
            usz row = to_matrix_id(i, j, mesh_size);
            matrix->data[CSC_idx] = 4.0;
            matrix->row_index[CSC_idx] = row;
            CSC_idx ++;

            //Close right
            if(i < mesh_size -1)
            {
                usz row = to_matrix_id(i+1,j, mesh_size);
                matrix->data[CSC_idx] = -1.0;
                matrix->row_index[CSC_idx] = row;
                CSC_idx ++;
            }
            //Far right
            if(j < mesh_size -1)
            {
                usz row = to_matrix_id(i,j+1, mesh_size);
                matrix->data[CSC_idx] = -1.0;
                matrix->row_index[CSC_idx] = row;
                CSC_idx ++;
            }

            matrix->col_index[col_idx] = CSC_idx;
            col_idx ++;
        }
    }
}
*/
void print_CSC(csc_matrix_t *matrix)
{
    assert(matrix->data != NULL);
 
    usz idx = 0;
    
    printf("matrix size : %ld\n", matrix->size);
    for(usz i = 0; i < matrix->col_index[matrix->size]; i++)
    {
        printf("%f ", matrix->data[i]);
    }
    printf("\n");
    printf("\n");

    for(usz i = 0; i < matrix->col_index[matrix->size]; i++)
    {
        printf("%d ", matrix->row_index[i]);
    }
    printf("\n");
    printf("\n");

    for(usz i = 0; i < matrix->size; i++)
    {
        printf("%d ", matrix->col_index[i]);
    }
    printf("\n");
    
    for(usz i = 0; i < matrix->size; i++)
    {
        //usz elems_in_row = matrix->row_index[i+1] - matrix->row_index[i]; 
        
        for(usz j = 0; j < matrix->size; j++)
        {
            //if(matrix->row_index[idx] == i && matrix->col_index[j] > 0)
            
            if (idx < matrix->size && i == matrix->row_index[idx] && matrix->col_index[j] > 0)
            {
                printf("%-2.3f ", matrix->data[idx]);
                //elems_in_row --;
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

void free_CSC(csc_matrix_t *matrix)
{
    assert(matrix != NULL);

    free(matrix->col_index);
    free(matrix->row_index);
    free(matrix->data);
}
