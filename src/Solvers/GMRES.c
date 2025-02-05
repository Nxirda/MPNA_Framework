#include "GMRES.h"

void arnoldi(matrix_t *matrix, matrix_t *Q, u64 k)
{
    vector_t q;
    vector_t h;

    q = A * Q(:,k);   // Krylov Vector
    
    for (usz i = 0; i < k; i++)       // Modified Gram-Schmidt, keeping the Hessenberg matrix
    {
        h[i] = q * Q(:,i);
        q = q - h[i] * Q(:,i);
    }
    h->data[k+1] = norm(q);
    q = q / h->data[k+1];
}

void GMRES_general(matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol)
{

    const usz N = x->size;
    f64 (*A)[N] = make_2D_span(f64, , matrix->data, N);        

    vector_t residual;
    allocate_vector(&residual, N);
    compute_residual_general(matrix, x, b, &residual);
 
    u64 norm_b = ;
    u64 norm_r = ;

    u64 error = norm_r / norm_b;

    vector_t sn;
    allocate_vector(&sn, max_iterations);
    
    vector_t cs;
    allocate_vector(&cs, max_iterations);

    vector_t e1;
    allocate_vector(&e1, max_iterations +1);
    e1->data[0] = 1;

    vector_t e;
    allocate_vector(&e, max_iterations);
    e->data[0] = error;

    vector_t beta;
    allocate_vector(&beta, max_iterations +1);
    beta = norm_r * e1;

    while((k < max_iterations) && (error < tol))
    {
        arnoldi_iteration
        apply_rotations

        // Update residual
        beta[k+1]   = -sn[k] * beta[k];
        beta[k]     = cs[k] * beta[k];
        error       = abs(beta[k +1]) / norm_b; 
        // Update error
        error[k+1] = error; 

        k++;
    }
    // Compute res
    //y = H \ beta;
    //x = x + Q * y;
}

void GMRES_csr(csr_matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol)
{
    return;
}
