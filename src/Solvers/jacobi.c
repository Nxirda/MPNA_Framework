#include "jacobi.h"

#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

void jacobi_general(matrix_t const *matrix, vector_t *vx, vector_t const *vb, u64 max_iterations)
{
    assert(matrix && vx && vb );

    const usz N = vx->size;

    vector_t swap_x;
    allocate_vector(&swap_x, N);
    
    f64 const *b = vb->data;

    f64 (*A)[N] = make_2D_span(f64, , matrix->values, N); 
    
    i8 converged = 0;
    u64 k = 0;

    while(!converged && k < max_iterations)
    {
        f64 *x_kp1  = swap_x.data;
        f64 *x      = vx->data;

        for(usz i = 0; i < N; i++)
        {
            f64 tmp_x = 0.0;
            f64 tmp_b = b[i];
            f64 tmp_a = 1.0/A[i][i];
           
            // Distinct loops instead of ifs : basic stuff
            for(usz j = 0; j < i; j++)
            {
                tmp_x += A[i][j] * x[j];
            }

            for(usz j = i+1; j < N;  j++)
            {
                tmp_x += A[i][j] * x[j];
            }

            x_kp1[i] = tmp_a * (tmp_b - tmp_x);    
        } 
        swap_vector(vx, &swap_x);
        k ++;
    }
    free_vector(&swap_x);
}

void jacobi_csr(csr_matrix_t const *matrix, vector_t *vx, vector_t const *vb, u64 max_iterations)
{
    const usz N = vx->size;
    vector_t swap_x;
    allocate_vector(&swap_x, N);

    f64 const *b = vb->data;
    u64 k = 0;
    i8 converged = 0;

    while(!converged && k < max_iterations)
    {
        usz idx = 0;
        f64 *x_kp1  = swap_x.data;
        f64 *x      = vx->data;

        for(u64 i = 0; i < N; i++)
        {
            usz row_elems = matrix->row_index[i+1] - matrix->row_index[i];
            f64 tmp_x = 0.0;
            f64 tmp_b = b[i];
            f64 tmp_a = 0.0; 

            for(usz j = 0; j < row_elems; j++)
            {
                usz curr_col = matrix->col_index[idx];
                
                if(curr_col != i)
                {
                    f64 Aij = matrix->values[idx];
                    tmp_x += Aij * x[curr_col];
                }
                else
                {
                    tmp_a = 1.0/matrix->values[idx];
                }
                idx ++;
            }
            x_kp1[i] = tmp_a * (tmp_b - tmp_x);
        }
        swap_vector(vx, &swap_x);
        k++;
    }
    free_vector(&swap_x);
}
