#include "jacobi.h"
#include <math.h>
#include <assert.h>
void jacobi_general(matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    assert(matrix && x && b );

    const usz N = x->size;

    vector_t swap_x;
    allocate_vector(&swap_x, N);
    
    f64 (*A)[N] = make_2D_span(f64, , matrix->data, N); 
    
    u8 converged = 0;
    u64 k = 0;
    f64 tmp_a = 0.0;
    f64 tmp_b = 0.0;
    f64 tmp_x = 0.0;
    f64 norm  = 0.0;
    while(!converged && k < max_iterations)
    {
        for(usz i = 0; i < N; i++)
        {
            tmp_x = 0.0;
            tmp_b = b->data[i];
            tmp_a = 1.0/A[i][i];
           
            // Distinct loops instead of ifs : basic stuff
            for(usz j = 0; j < i; j++)
            {
                tmp_x += A[i][j] * x->data[j];
            }

            for(usz j = i+1; j < N;  j++)
            {
                tmp_x += A[i][j] * x->data[j];
            }

            swap_x.data[i] = tmp_a * (tmp_b - tmp_x);    
        } 
        norm = 0.0;
        for(usz i = 0; i < N; i++)
        {
            norm += (swap_x.data[i] - x->data[i]) * (swap_x.data[i] - x->data[i]);
        }
        norm = sqrt(norm);

        if(norm < tol)
        {
            converged = 1;
        }
        swap_vector(x, &swap_x);
        k ++;
    }
    free_vector(&swap_x);
}

void jacobi_csr(csr_matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    const usz N = x->size;
    i8 converged = 0;

    vector_t swap_x;
    allocate_vector(&swap_x, N);
    
    u64 k           = 0;
    usz idx         = 0;
    usz nb_elems    = 0;
    usz curr_col    = 0;

    f64 tmp_a   = 0.0;
    f64 tmp_b   = 0.0;
    f64 tmp_x   = 0.0;
    f64 Aij     = 0.0;
    f64 norm    = 0.0;
    while(!converged && k < max_iterations)
    {
        idx = 0;

        for(usz i = 0; i < N; i++)
        {
            tmp_a = 0.0;
            tmp_b = b->data[i];
            tmp_x = 0.0; 
            nb_elems = matrix->row_index[i+1] - matrix->row_index[i];
            
            for(usz j = 0; j < nb_elems; j++, idx++)
            {
                curr_col = matrix->col_index[idx];
                Aij = matrix->data[idx];
                
                if(curr_col != i)
                {
                    tmp_x += Aij * x->data[curr_col];
                }
                else
                {
                    tmp_a = 1.0/Aij;
                }
            }
            swap_x.data[i] = tmp_a * (tmp_b - tmp_x);
        }

        norm = 0.0;
        for(usz i = 0; i < N; i++)
        {
            norm += (swap_x.data[i] - x->data[i]) * (swap_x.data[i] - x->data[i]);
        }
        norm = sqrt(norm);

        if(norm < tol)
        {
            converged = 1;
        }
        //converged = allclose_vector(x, &swap_x, 1e-08, 1e-05);
        swap_vector(x, &swap_x);
        k++;
    }
    free_vector(&swap_x);
}
