#include "GMRES.h"
#include "linear_algebra.h"
#include <lapacke.h>
#include <math.h>

/*
 * S/o to this guy that is an absolute unit, basically every thing you need
 * he did it
 *
 * https://people.math.sc.edu/Burkardt/c_src/qr_solve/qr_solve.c
 *
 * John Burkardt's the goat
 *
 * Now the problem lies in the fact i did not understood anything, so after beeing
 * defeated by this approach I decided to use the wikipedia much simpler approach
 * provided here : 
 *
 * https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
 *
 * */

// Concretely solves Ax = b into r
void compute_residual_general_2(matrix_t const *matrix, vector_t const *x, vector_t const *b, vector_t *r)
{
    const usz N = x->size;
    f64 (*A)[matrix->dim_y] = make_2D_span(f64, , matrix->data, matrix->dim_y);

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
}

void arnoldi_step(matrix_t const *mat_A, matrix_t *mat_Q, matrix_t *mat_H, u64 k)
{
    f64 (*A)[mat_A->dim_y]  = make_2D_span(f64, , mat_A->data, mat_A->dim_y);
    f64 (*H)[mat_H->dim_y]  = make_2D_span(f64, , mat_H->data, mat_H->dim_y);
    f64 (*Q)[mat_Q->dim_y]  = make_2D_span(f64, , mat_Q->data, mat_Q->dim_y);

    // qk = A * qk-1 | Krylov-Vector
    for(usz i = 0; i < mat_Q->dim_x; i++)
    {
        Q[i][k] = 0.0;
        for(usz j = 0; j < mat_A->dim_y; j++)
        {
            Q[i][k] += A[i][j] * Q[j][k-1]; 
        }
    }
   
    // Modified Gram Schmidt keeping the hessenberg matrix
    for(usz j = 0; j < k; j++)
    {
        // Hj,k-1 = qj * qk
        H[j][k-1] = 0.0;
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            H[j][k-1] += Q[i][j] * Q[i][k];
        }

        // qk = qk - Hj,k-1 * qj
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            Q[i][k] -= H[j][k-1] * Q[i][j];
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
    if(norm_q != 0.0)
    {
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            Q[i][k] *= (1.0/norm_q);
        }
    }
}

void rotations(usz k, f64 a, f64 b, vector_t *cs, vector_t *sn)
{
    f64 a_2 = a * a;
    f64 b_2 = b * b;
    f64 t = sqrt(a_2 + b_2);
    cs->data[k] = a / t;
    sn->data[k]=  b / t;
}

void apply_rotations(usz k, matrix_t *mat_H, vector_t *cs, vector_t *sn)
{
    //usz curr_k = k-1;
    f64(* H)[mat_H->dim_y] = make_2D_span(f64, , mat_H->data, mat_H->dim_y);
    for(usz i = 0; i < k-1; i++)
    {
        f64 tmp     =  cs->data[i] * H[i][k-1] + sn->data[i] * H[i+1][k-1];
        H[i+1][k-1] = -sn->data[i] * H[i][k-1] + cs->data[i] * H[i+1][k-1];
        H[i][k-1]   = tmp;
    }

    rotations(k-1, H[k-1][k-1], H[k][k-1], cs, sn);

    H[k-1][k-1] = cs->data[k-1] * H[k-1][k-1] + sn->data[k-1] * H[k][k-1];
    H[k][k-1] = 0.0;
}

// Solve least-squares problem ||g - H * y ||, using back substitution 
void solve_ls(usz k, vector_t *g, matrix_t *mat_H, vector_t *y)
{
    f64(* H)[mat_H->dim_y] = make_2D_span(f64, ,mat_H->data, mat_H->dim_y);
    for(usz i = k-1; i > 0; i--)
    {
        f64 sum = g->data[i];
        for(usz j = i+1; j < k; j++)
        {
            sum -= H[i][j] * y->data[j];
        }
        y->data[i] = sum * (1.0/H[i][i]);
    }
}

usz GMRES_general(matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol)
{
    usz N = b->size;
    const u64 size = max_iterations; // Concretely this is the maximum Krylov subspace dim

    matrix_t Q;
    allocate_matrix(N, size+1, &Q);
    fill_matrix(&Q, 0.0);

    matrix_t H;
    allocate_matrix(size+1, size, &H);
    fill_matrix(&H, 0.0);

    vector_t residual;
    allocate_vector(&residual, N);
    compute_residual_general_2(matrix, x, b, &residual);
    
    f64 norm = dot_product(&residual, &residual);
    f64 beta = sqrt(norm);

    f64 (*Q_span)[Q.dim_y] = make_2D_span(f64, ,Q.data, Q.dim_y);
    f64 (*H_span)[H.dim_y] = make_2D_span(f64, ,H.data, H.dim_y);

    for(usz i = 0; i < Q.dim_x; i++)
        Q_span[i][0] = residual.data[i] / beta;

    /*vector_t y;
    allocate_vector(&y, size+1);
    init_constant_vector(&y, 0.0);*/
    
    vector_t tmp_x;
    allocate_vector(&tmp_x, x->size);
    //vector_t cs;
    //vector_t sn;

    /*allocate_vector(&cs, size);
    init_constant_vector(&cs, 0.0);

    allocate_vector(&sn, size);
    init_constant_vector(&sn, 0.0);*/

    vector_t g;
    allocate_vector(&g, size+1);
    init_constant_vector(&g, 0.0);
    g.data[0] = beta;

    f64 *H_1d = (f64 *)malloc((size+1) * max_iterations * sizeof(f64));
    
    usz k = 1;
    while(k < max_iterations)
    {
        arnoldi_step(matrix, &Q, &H, k);
       
        init_constant_vector(&g, 0.0);
        g.data[0] = beta;


        for (int i = 0; i < size+1; ++i) {
            for (int j = 0; j < k; ++j) {
                H_1d[j*(size+1) + i] = H_span[i][j];
            }
        }

        LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', size+1, k, 1, H_1d, size+1, g.data, size+1);
    
        //copy_vector(&g, &y);
        //print_vector(&g);
        
        for(usz i = 0; i < Q.dim_x; i++)
        {
            tmp_x.data[i] = x->data[i];
            for(usz j = 0; j < Q.dim_y -1; j++)
            {
                tmp_x.data[i] += Q_span[i][j] * g.data[j];
            }
        }

        compute_residual_general_2(matrix, &tmp_x, b, &residual);
        norm = dot_product(&residual, &residual);
        norm = sqrt(norm);
        
        //apply_rotations(k, &H, &cs, &sn);
        //g.data[k] = -sn.data[k-1] * g.data[k-1];
        //g.data[k-1] = cs.data[k-1] * g.data[k-1];
        // Might lack a condition here tbh
        if(norm < tol)
        {
            break;
        }

        k++;
    }
   
    // g - H*y;
    //solve_ls(k, &g, &H, &y); 
    // x = x + Q * y || consider only Q until K
    for(usz i = 0; i < Q.dim_x; i ++)
    {
        f64 sum = 0.0;
        for(usz j = 0; j < k; j++)
        {
            sum += Q_span[i][j] * g.data[j]; 
        }
        x->data[i] = sum;
    }

    free(H_1d);
    free_matrix(Q);
    free_matrix(H);
    free_vector(&residual);
    //free_vector(&y);
    free_vector(&tmp_x);
    free_vector(&g);
    return k;
}

void arnoldi_step_csr(csr_matrix_t const *A, matrix_t *mat_Q, matrix_t *mat_H, u64 k)
{
    //f64 (*A)[mat_A->dim_y]  = make_2D_span(f64, , mat_A->data, mat_A->dim_y);
    f64 (*H)[mat_H->dim_y]  = make_2D_span(f64, , mat_H->data, mat_H->dim_y);
    f64 (*Q)[mat_Q->dim_y]  = make_2D_span(f64, , mat_Q->data, mat_Q->dim_y);

    // qk = A * qk-1 | Krylov-Vector
    for(usz i = 0; i < mat_Q->dim_x; i++)
    {
        Q[i][k] = 0.0;
        //for(usz j = 0; j < mat_A->dim_y; j++)
        for(usz j = A->row_index[i]; j < A->row_index[i+1]; j++)
        {
            //Q[i][k] += A[i][j] * Q[j][k-1]; 
            Q[i][k] += A->data[j] * Q[A->col_index[j]][k-1];
        }
    }
   
    // Modified Gram Schmidt keeping the hessenberg matrix
    for(usz j = 0; j < k; j++)
    {
        // Hj,k-1 = qj * qk
        H[j][k-1] = 0.0;
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            H[j][k-1] += Q[i][j] * Q[i][k];
        }

        // qk = qk - Hj,k-1 * qj
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            Q[i][k] -= H[j][k-1] * Q[i][j];
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
    if(norm_q != 0.0)
    {
        for(usz i = 0; i < mat_Q->dim_x; i++)
        {
            Q[i][k] *= (1.0/norm_q);
        }
    }
}
//
void compute_residual_csr_2(const csr_matrix_t *matrix, const vector_t *x, const vector_t *b, vector_t *r)
{
    const usz N = x->size;
    
    for(usz i = 0; i < N; i++)
    {
        f64 r_i = 0.0; 
        f64 b_i = b->data[i];

        for(usz j = matrix->row_index[i]; j < matrix->row_index[i+1]; j++)
        {
            r_i += matrix->data[j] * x->data[matrix->col_index[j]];
        }
        
        r->data[i] = b_i - r_i;
    }
}

usz GMRES_csr(csr_matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol)
{
    usz N = b->size;
    const u64 size = max_iterations; // Concretely this is the maximum Krylov subspace dim

    matrix_t Q;
    allocate_matrix(N, size+1, &Q);
    fill_matrix(&Q, 0.0);

    matrix_t H;
    allocate_matrix(size+1, size, &H);
    fill_matrix(&H, 0.0);

    vector_t residual;
    allocate_vector(&residual, N);
    compute_residual_csr_2(matrix, x, b, &residual);
    
    f64 norm = dot_product(&residual, &residual);
    f64 beta = sqrt(norm);

    f64 (*Q_span)[Q.dim_y] = make_2D_span(f64, ,Q.data, Q.dim_y);
    f64 (*H_span)[H.dim_y] = make_2D_span(f64, ,H.data, H.dim_y);

    for(usz i = 0; i < Q.dim_x; i++)
        Q_span[i][0] = residual.data[i] / beta;

    /*vector_t y;
    allocate_vector(&y, size+1);
    init_constant_vector(&y, 0.0);*/
    
    vector_t tmp_x;
    allocate_vector(&tmp_x, x->size);
    
    vector_t g;
    allocate_vector(&g, size+1);
    init_constant_vector(&g, 0.0);
    g.data[0] = beta;

    f64 *H_1d = (f64 *)malloc((size+1)* max_iterations * sizeof(f64));
    
    usz k = 1;
    while(k < max_iterations)
    {
        arnoldi_step_csr(matrix, &Q, &H, k);
       
        init_constant_vector(&g, 0.0);
        g.data[0] = beta;


        for (int i = 0; i < size+1; ++i) {
            for (int j = 0; j < k; ++j) {
                H_1d[j*(size+1) + i] = H_span[i][j];
            }
        }

        LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', size+1, k, 1, H_1d, size+1, g.data, size+1);
        
        for(usz i = 0; i < Q.dim_x; i++)
        {
            tmp_x.data[i] = x->data[i];
            for(usz j = 0; j < Q.dim_y -1; j++)
            {
                tmp_x.data[i] += Q_span[i][j] * g.data[j];
            }
        }

        compute_residual_csr_2(matrix, &tmp_x, b, &residual);
        norm = dot_product(&residual, &residual);
        norm = sqrt(norm);
        
        if(norm < tol)
        {
            break;
        }

        k++;
    }
   
    for(usz i = 0; i < Q.dim_x; i ++)
    {
        f64 sum = 0.0;
        for(usz j = 0; j < k; j++)
        {
            sum += Q_span[i][j] * g.data[j]; 
        }
        x->data[i] += sum;
    }

    free(H_1d);
    free_matrix(Q);
    free_matrix(H);
    free_vector(&g);
    free_vector(&residual);
    free_vector(&tmp_x);
    return k;
}
