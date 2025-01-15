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
void gauss_seidel_general(matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    assert(matrix && x && b );

    const usz N = x->size;
    f64 (*A)[N] = make_2D_span(f64, , matrix->data, N); 
    
    i8 converged = 0;
    u64 k = 0; 
    f64 tmp_x = 0.0;
    f64 tmp_a = 0.0;
    f64 tmp_b = 0.0;
    f64 norm  = 0.0;
    f64 new   = 0.0;
    while(!converged && k < max_iterations)
    {
        norm = 0.0;
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
            new = tmp_a * (tmp_b - tmp_x); 
            norm += (new - x->data[i]) * (new - x->data[i]);
            x->data[i] = new;    
        }

        norm = sqrt(norm);
        if(norm < tol)
        {
            converged = 1;
        }
        k ++;
    }
}

void gauss_seidel_csr(csr_matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations, f64 tol)
{
    const usz N = x->size; 
    i8 converged = 0;
    
    u64 k           = 0;
    usz idx         = 0;
    usz nb_elems    = 0;
    usz curr_col    = 0;

    f64 tmp_a   = 0.0;
    f64 tmp_b   = 0.0;
    f64 tmp_x   = 0.0;
    f64 Aij     = 0.0;
    f64 new     = 0.0;
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
                else // diagonal elem
                {
                    tmp_a = 1.0/Aij;
                }
            }
            new = tmp_a * (tmp_b - tmp_x); 
            norm += (new - x->data[i]) * (new - x->data[i]);
            x->data[i] = new;    
        }
        norm = sqrt(norm);
        if(norm < tol)
        {
            converged = 1;
        }
        k++;
    }
}
