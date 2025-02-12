#include "GMRES.h"

#include <math.h>
/*
void compute_residual_general(matrix_t const *matrix, vector_t const *x, vector_t const *b, vector_t *r)
{
    const usz N = x->size;
    f64 (*A)[N] = make_2D_span(f64, , matrix->data, N);

    for(usz i = 0; i < N; i++)
    {
        f64 r_i = 0.0;
        f64 b_i = b->data[i];

        for(usz j = 0; j < N; j++)
        {
            r_i += A[i][j] * x->data[j]; 
        }
        
        r->data[i] = b_i - r_i;
    }
}*/

inline void arnoldi(matrix_t const *mat_A, matrix_t *mat_Q, matrix_t *mat_H, u64 lim_k)
{
    f64 (*A)[mat_A->size]   = make_2D_span(f64, , mat_A->data, mat_A->size);
    // Would help tbh if those were matrices of vectors
    f64 (*H)[mat_H->dim_y]  = make_2D_span(f64, , mat_H->data, mat_H->dim_y);
    f64 (*Q)[mat_Q->dim_y]  = make_2D_span(f64, , mat_Q->data, mat_Q->dim_y);

    f64 tmp_q = 0.0;
    for (usz k = 1; k < lim_k; k++) // Modified Gram-Schmidt, keeping the Hessenberg matrix
    {
        // qk = A * qk-1
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            Q[i][k] = 0;
            for(usz j = 0; j < mat_A->dim_y; j++)
            {
                Q[i][k] += A[i][j] * Q[j][k-1]; 
            }
        }
        
        for(usz j = 0; j < k; j++)
        {
            // Hj,k-1 = qj * qk
            for(usz i = 0; i < mat_Q->dim_x; i++ )
            {
                H[j][k-1] = Q[i][j] * Q[i][k];
            }

            // qk = qk - Hj,k-1 * qj
            for(usz i = 0; i < mat_Q->dim_x; i++)
            {
                Q[i][k] = Q[i][k] - H[j][k-1] * Q[i][j];
            }
        }

        f64 norm_q = 0.0;
        
        // Hk,k-1 = norm(qk)
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            norm_q += Q[i][k] * Q[i][k];
        }

        norm_q = sqrt(norm_q);
        H[k][k-1] = norm_q;
        
        // qk = qk/ Hk,k-1
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            Q[i][k] = Q[i][k] * (1.0/norm_q);
        }
    }
}

void GMRES_general(matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol)
{
    const usz N = x->size;
    //f64 (*A)[N] = make_2D_span(f64, , matrix->data, N);        

    matrix_t Q;
    allocate_rectangle_matrix(b->size, N+1, Q);
    f64 (*Q_span)[N+1] = make_2D_span(f64, ,Q->data, N+1);


    matrix_t H;
    allocate_rectangle_matrix(N+1, N, H);
    fill_matrix(0, H.size);

    vector_t residual;
    allocate_vector(&residual, N);
    compute_residual_general(matrix, x, b, &residual);
    
    u64 norm_r = dot_product(&residual, &residual);
    
    for(usz j = 0; j < Q.dim_x; j++)
    {
        Q_span[j][0] = residual.data[j] / norm_r;
    }

    usz k = 0;
    while((k < max_iterations))// && (error < tol))
    {
        arnoldi(matrix, &Q, &H, lim_k);
        // Update residual
        k++;
    }
}

/*void GMRES_csr(csr_matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol)
{
    return;
}*/
