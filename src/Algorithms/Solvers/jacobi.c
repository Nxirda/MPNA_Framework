#include "jacobi.h"
#include <math.h>
#include <assert.h>

usz jacobi_general(matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    assert(matrix && x && b );

    const usz N = x->size;

    vector_t swap_x;
    allocate_vector(&swap_x, N);
    init_constant_vector(&swap_x, 0.0);
    
    f64 (*A)[matrix->dim_y] = make_2D_span(f64, , matrix->data, matrix->dim_y); 
    
    u64 k       = 0;
    f64 diag    = 0.0;
    f64 sum     = 0.0;
    f64 norm    = 0.0;
    while(k < max_iterations)
    {
        for(usz i = 0; i < N; i++)
        {
            sum = 0.0;
            diag = A[i][i];
           
            // Distinct loops instead of ifs : basic stuff
            for(usz j = 0; j < i; j++)
            {
                sum += A[i][j] * x->data[j];
            }

            for(usz j = i+1; j < N;  j++)
            {
                sum += A[i][j] * x->data[j];
            }

            swap_x.data[i] = (b->data[i] - sum) * (1.0/diag);    
        } 
        norm = 0.0;
        for(usz i = 0; i < N; i++)
        {
            norm += (swap_x.data[i] - x->data[i]) * (swap_x.data[i] - x->data[i]);
        }
        norm = sqrt(norm);

        if(norm < tol)
        {
            free_vector(&swap_x);
            return k;
        }
        swap_vector(x, &swap_x);
        k ++;
    }
    free_vector(&swap_x);
    return max_iterations;
}

usz jacobi_csr(csr_matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    const usz N = x->size;

    vector_t swap_x;
    allocate_vector(&swap_x, N);
    
    usz k           = 0;
    usz curr_col    = 0;

    f64 diag    = 0.0;
    f64 sum     = 0.0;
    f64 norm    = 0.0;
    while(k < max_iterations)
    {
        for(usz i = 0; i < N; i++)
        {
            diag = 0.0;
            sum = 0.0; 
            
            for(usz j = matrix->row_index[i]; j < matrix->row_index[i+1]; j++)
            {
                curr_col = matrix->col_index[j];
                
                if(curr_col != i)
                {
                    sum += matrix->data[j] * x->data[curr_col];
                }
                else
                {
                    diag = matrix->data[j];
                }
            }
            swap_x.data[i] = (b->data[i] - sum) * (1.0/diag);
        }

        norm = 0.0;
        for(usz i = 0; i < N; i++)
        {
            norm += (swap_x.data[i] - x->data[i]) * (swap_x.data[i] - x->data[i]);
        }
        norm = sqrt(norm);

        if(norm < tol)
        {
            free_vector(&swap_x);
            return k;
        }
        swap_vector(x, &swap_x);
        k++;
    }
    free_vector(&swap_x);
    return max_iterations;
}
