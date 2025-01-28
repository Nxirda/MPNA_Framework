#include "conjugate_gradient.h"

//
void compute_residual_general(matrix_t const *matrix, vector_t const *x, vector_t const *b, vector_t *r)
{
    const usz N = x->size;
    f64 (*A)[N] = make_2D_span(f64, , matrix->data, N);

    for(usz i = 0; i < N; i++)
    {
        f64 x_i = x->data[i];
        f64 r_i = 0.0;
        f64 b_i = b->data[i];
        for(usz j = 0; j < N; j++)
        {
            r_i += A[i][j] * x_i;
        }
        
        r->data[i] = b_i - r_i;
    }
}

//
void conjugate_gradient_general(matrix_t const *matrix, vector_t *x, 
                vector_t const *b, u64 max_iterations, f64 tol)
{
    u64 k = 0;
    const usz N = x->size;
    f64 (*A)[N] = make_2D_span(f64, , matrix->data, N);

    vector_t residual;
    allocate_vector(&residual, N);
    compute_residual_general(matrix, x, b, &residual);
    
    vector_t direction;
    allocate_vector(&direction, N);
    copy_vector(&residual, &direction);

    vector_t q;
    allocate_vector(&q, N);

    f64 d_new = dot_product(&residual, &residual);
    f64 d_old = d_new;
    
    f64 tmp_q = 0.0;
    f64 tmp_d = 0.0;
    f64 alpha = 0.0;
    f64 beta  = 0.0;

    while((k < max_iterations) && (d_new > (tol * tol) * d_old))
    {
        for(usz i = 0; i < N; i++)
        {
            tmp_q   = 0.0;
            tmp_d     = direction.data[i];    
        
            for(usz j = 0; j < N; j++)
            {
                tmp_q += A[i][j] * tmp_d;      
            }

            q.data[i] = tmp_q;
        }
       
        f64 tmp = dot_product(&direction, &q);
        alpha = d_new / tmp;
        
        daxpy(alpha, &direction, x, x);
        //x->data[i] += alpha * tmp_d;

        //residual.data[i] += -(alpha * tmp_q);
        if((k%50) == 0)
            compute_residual_general(matrix, x, b, &residual);
        else
            daxpy(-alpha, &q, &residual, &residual);

        d_old = d_new;
        d_new = dot_product(&residual, &residual);
        
        beta = d_new / d_old; 
        // residual = beta * direction + direction 
        daxpy(beta, &direction, &residual, &direction);
        k ++;
    }
}

//
void compute_residual_csr(const csr_matrix_t *matrix, const vector_t *x, const vector_t *b, vector_t *r)
{
    const usz N = x->size;
    
    usz idx         = 0;
    usz nb_elems    = 0;

    for(usz i = 0; i < N; i++)
    {
        f64 x_i = x->data[i];
        f64 r_i = b->data[i];

        nb_elems = matrix->row_index[i+1] - matrix->row_index[i];

        for(usz j = 0; j < nb_elems; j++, idx++)
        {
            r_i -= matrix->data[idx] * x_i;
        }
        
        r->data[i] = r_i;
    }
}

//
void conjugate_gradient_csr(csr_matrix_t const *matrix, vector_t *x, 
            vector_t const *b, u64 max_iterations, f64 tol)
{
    u64 k = 0;
    const usz N = x->size;

    vector_t residual;
    allocate_vector(&residual, N);
    compute_residual_csr(matrix, x, b, &residual);
    
    vector_t direction;
    allocate_vector(&direction, N);
    copy_vector(&residual, &direction);

    f64 d_new = dot_product(&residual, &residual);
    f64 d_old = d_new;

    usz idx         = 0;
    usz nb_elems    = 0;
    
    f64 tmp_q = 0.0;
    f64 tmp_d = 0.0;
    f64 alpha = 0.0;
    f64 beta  = 0.0;
   
    while(k < max_iterations && d_new > (tol * tol) * d_old)
    {
        idx = 0;
        for(usz i = 0; i < N; i++)
        {
            tmp_q = 0.0;
            tmp_d = direction.data[i]; 

            for(usz j = matrix->row_index[i]; j < matrix->row_index[i+1] ; j++)
            {
                tmp_q += matrix->data[j] * tmp_d;      
            }

            alpha = d_new / (tmp_d* tmp_q);
            
            x->data[i] += alpha * tmp_d;
        }
        
        compute_residual_csr(matrix, x, b, &residual);
        
        d_old = d_new;
        d_new = dot_product(&residual, &residual);
        
        beta = d_new / d_old; 
        
        // residual = residual + beta * direction 
        daxpy(beta, &direction, &residual, &direction);

        k ++;
    }
}
