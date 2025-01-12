#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "general.h"

/*
    Given a mesh size n, allocate the corresponding 
    matrix of size N * N
    where N = n^2
*/
void allocate_matrix(usz mesh_size, matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);
    usz N = mesh_size * mesh_size;

    matrix->size = N;
    matrix->data = (f64 *)calloc(N * N, sizeof(f64));
    //printf("Allocated size : %ld\n", N * N * sizeof(f64));
}

void fill_matrix(usz mesh_size, matrix_t matrix)
{
    assert(matrix.data != NULL && mesh_size > 0);

    const usz N = matrix.size;
    f64 (* tmp_matrix)[N] = make_2D_span(f64, ,matrix.data, N);
    
    // This order is due to how we parse the mesh into a matrix
    for(usz j = 0; j < mesh_size; j++)
    {
        for(usz i = 0; i < mesh_size; i++)
        {
            // Corresponding coord in matrix
            const usz k = to_matrix_id(i,j,mesh_size);

            //Far left 
            if(j > 0)
                tmp_matrix[k][to_matrix_id(i,j-1,mesh_size)] = -1.0;
            
            //Close left
            if(i > 0)
                tmp_matrix[k][to_matrix_id(i-1,j,mesh_size)] = -1.0;

            tmp_matrix[k][k] = 4.0;

            //Close right
            if(i < mesh_size -1)
                tmp_matrix[k][to_matrix_id(i+1,j,mesh_size)] = -1.0;
             
            //Far right
            if(j < mesh_size -1)
                tmp_matrix[k][to_matrix_id(i,j+1,mesh_size)] = -1.0;
        }
    }
}

void print_matrix(matrix_t matrix)
{
    usz N = matrix.size;
    f64 (*tmp_matrix)[N] = make_2D_span(f64, , matrix.data, N);

    for(usz i = 0; i < N; i++)
    {
        for(usz j = 0; j < N; j++)
        {
            if(tmp_matrix[i][j] == 0)
                printf("__.___ ");
            else
                printf("%-2.3f ", tmp_matrix[i][j]);
        }
        printf("\n");
    }
}

void free_matrix(matrix_t matrix)
{
    assert(matrix.data != NULL);
    free(matrix.data);
}
