#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "general.h"

/*
    Given a mesh size n, allocate the corresponding 
    matrix of size N * N
    where N = n^2
*/
/*void allocate_matrix(usz mesh_size, matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);
    usz N = mesh_size * mesh_size;
    
    matrix->size = N;
    matrix->data = (f64 *)calloc(N * N, sizeof(f64));
    //printf("Allocated size : %ld\n", N * N * sizeof(f64));
}*/

void allocate_matrix(usz size_x, usz size_y, matrix_t *matrix)
{
    assert(matrix != NULL && size_x > 0 && size_y > 0);
    const usz N = size_x * size_y;

    matrix->size = N;
    matrix->dim_x = size_x;
    matrix->dim_y = size_y;
    matrix->data = (f64 *)malloc(N * sizeof(f64));
}

void fill_matrix(matrix_t *matrix, f64 value)
{
    assert(matrix->data != NULL);
    
    f64(* span)[matrix->dim_y] = make_2D_span(f64, ,matrix->data, matrix->dim_y);

    for(usz i = 0; i < matrix->dim_x; i++)
    {
        for(usz j = 0; j < matrix->dim_y; j++)
        {
            span[i][j] = value;
        }
    }
}

/*void fill_matrix(usz mesh_size, matrix_t matrix)
{
    assert(matrix.data != NULL && mesh_size > 0);

    const usz N = matrix.size;
    f64 (* span)[N] = make_2D_span(f64, ,matrix.data, N);
    
    // This order is due to how we parse the mesh into a matrix
    for(usz j = 0; j < mesh_size; j++)
    {
        for(usz i = 0; i < mesh_size; i++)
        {
            // Corresponding coord in matrix
            const usz k = to_matrix_id(i,j,mesh_size);

            //Far left 
            if(j > 0)
                span[k][to_matrix_id(i,j-1,mesh_size)] = -1.0;
            
            //Close left
            if(i > 0)
                span[k][to_matrix_id(i-1,j,mesh_size)] = -1.0;

            span[k][k] = 4.0;

            //Close right
            if(i < mesh_size -1)
                span[k][to_matrix_id(i+1,j,mesh_size)] = -1.0;
             
            //Far right
            if(j < mesh_size -1)
                span[k][to_matrix_id(i,j+1,mesh_size)] = -1.0;
        }
    }
}
*/

void print_matrix(matrix_t matrix)
{
    //usz N = matrix.size;
    f64 (*span)[matrix.dim_y] = make_2D_span(f64, , matrix.data, matrix.dim_y);

    for(usz i = 0; i < matrix.dim_x; i++)
    {
        for(usz j = 0; j < matrix.dim_y; j++)
        {
            if(span[i][j] == 0)
                printf("__.___ ");
            else
                printf("%-2.3f ", span[i][j]);
        }
        printf("\n");
    }
}

void free_matrix(matrix_t matrix)
{
    assert(matrix.data != NULL);
    free(matrix.data);
}
