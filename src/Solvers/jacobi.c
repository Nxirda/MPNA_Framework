#include "jacobi.h"

#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

void jacobi_general(matrix_t const *matrix, vector_t *vx, vector_t const *vb, u64 max_iterations)
{
    assert(matrix != NULL && vx != NULL && vb != NULL);

    const usz N = matrix->size;

    vector_t swap_x;
    allocate_vector(&swap_x, N);
    
    f64 *x_kp1 = swap_x.data;
    f64 *x = vx->data;
    f64 *b = vb->data;

    f64 (*A)[N] = make_2D_span(f64, , matrix->values, N); 
    
    i8 converged = false;
    u64 k = 0;

    while(!converged && k < max_iterations)
    {
        for(usz i = 0; i < N; i++)
        {
            f64 tmp_x = 0.0;
            f64 tmp_b = b[i];
            f64 tmp_a = 1.0/A[i][i];
            
            for(usz j = 0; j < N; j++)
            {
                if(i == j) continue;
                tmp_x += A[i][j] * x[j];
            }
            x_kp1[i] = tmp_a * (tmp_b - tmp_x);    
        }
        k ++;
       
        f64 *tmp = swap_x.data;
        swap_x.data = vx->data;
        vx->data = tmp;
        
        // Temporary ofc
        converged = (k == 5 ? true : false);
    }

    free(swap_x.data);
}

void jacobi_csr(csr_matrix_t *matrix, f64 *x, f64 *b, usz size)
{
    return;
}
