#include "gauss_seidel.h"
#include <math.h>
#include <assert.h>

/*
    TO parallelize gauss seidel, we need to sweep the diagonal,
    elements below are done before and up top of the diagonal
    will be done in the next iter
    
    This is called a parallelisation by front

    Can we overlap the sweeps (sort of) ?
    Yes, we can do a train of weave to sweep the matrix
    globally this is a pipeline parallelism

    We can use omp tasks for this purpose but not efficient
    the scheduler will hurt really bad but it s a good way
    of understading and starting with omp tasks

    But we dont have "enough parallelism" -> red black gauss seidel
*/
usz gauss_seidel_general(matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    assert(matrix && x && b );

    const usz N = x->size;
    f64 (*A)[matrix->dim_y] = make_2D_span(f64, , matrix->data, matrix->dim_y); 
    
    u64 k = 0; 
    f64 sum = 0.0;
    f64 diag = 0.0;
    f64 norm  = 0.0;
    f64 new   = 0.0;
    while(k < max_iterations)
    {
        norm = 0.0;
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

            new = (b->data[i] - sum) * (1.0/diag); 
            norm += (new - x->data[i]) * (new - x->data[i]);
            x->data[i] = new;    
        }

        norm = sqrt(norm);
        if(norm < tol)
        {
            return k;
        }
        k ++;
    }
    return max_iterations;
}

usz gauss_seidel_csr(csr_matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    const usz N = x->size; 
    
    usz k           = 0;
    usz curr_col    = 0;

    f64 diag    = 0.0;
    f64 sum     = 0.0;
    f64 new     = 0.0;
    f64 norm    = 0.0;
    while(k < max_iterations)
    {
        norm = 0.0;
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
            new = (b->data[i] - sum) * (1.0/diag); 
            norm += (new - x->data[i]) * (new - x->data[i]);
            x->data[i] = new;    
        }
        norm = sqrt(norm);
        if(norm < tol)
        {
            return k;
        }
        k++;
    }
    return max_iterations;
}
