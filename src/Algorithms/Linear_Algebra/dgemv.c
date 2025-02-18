#include "dgemv.h" 

// r = alpha*A*x + beta*b
void general_mv(f64 alpha, matrix_t *A, vector_t *x, f64 beta, vector_t *b, vector_t *r)
{
    const usz N = x->size;
    f64 (*A_span)[A->dim_y] = make_2D_span(f64, , A->data, A->dim_y);

    for(usz i = 0; i < N; i++)
    {
        f64 r_i = 0.0;
        f64 b_i = b->data[i];

        for(usz j = 0; j < N; j++)
        {
            r_i += alpha * A_span[i][j] * x->data[j]; 
        }
        
        r->data[i] = beta * b_i + r_i;
    }
}

void csr_mv(f64 alpha, csr_matrix_t *A, vector_t *x, f64 beta, vector_t *b, vector_t *r)
{
    const usz N = x->size;

    for(usz i = 0; i < N; i++)
    {
        f64 r_i = 0.0;
        f64 b_i = b->data[i];

        for(usz j = A->row_index[i]; j < A->row_index[i+1]; j++)
        {
            r_i += alpha * A->data[j] * x->data[A->col_index[j]]; 
        }
        
        r->data[i] = beta * b_i + r_i;
    }    
}
