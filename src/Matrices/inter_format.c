#include "inter_format.h"

#include <stdlib.h>
#include <stdio.h>
usz NNZ_in_general(matrix_t const *a)
{
    const usz N = a->size;
    usz NNZ = 0;
    f64 (* a_span)[N] = make_2D_span(f64, , a->data, N);

    for(usz i = 0; i < N; i++)
    {
        for(usz j = 0; j < N; j++)
        {
            if(a_span[i][j] != 0.0)
            {
                NNZ +=1;
            }
        }
    }
    return NNZ;
}

// Meh is maybe an MPI Version tbh
void coo_to_csr(coo_matrix_t *matrix, csr_matrix_t *target)
{
    usz nnz = matrix->nnz;
    usz dim_x = 0, dim_y = 0;

    usz base_row_idx = matrix->row_index[0];
    usz base_col_idx = matrix->col_index[0];

    dim_x = matrix->dim_x; // We keep the same ammount of elemns per row
    dim_y = matrix->dim_y;
    
    allocate_CSR(dim_x, dim_y, nnz, target);
    
    for(usz i = 0; i < dim_y+1; i++)
        target->row_index[i] = 0;

    for(usz i = 0; i < nnz; i++)
    {
        target->data[i] = matrix->data[i];
        target->col_index[i] = matrix->col_index[i];
        target->row_index[matrix->row_index[i]+1]++;
    }

    for(usz i = 0; i < dim_y; i++)
        target->row_index[i+1] += target->row_index[i];
}
/*matrix_t csr_to_general(csr_matrix_t const *a)
{

    const usz N = a->size;

    u64 idx = 0;
    usz row_elems = 0;
    f64 *new_data = (f64 *)malloc(N * N * sizeof(f64));
    f64 (* a_span)[N] = make_2D_span(f64, , new_data, N);

    for(usz i = 0; i < N; i++)
    {
        row_elems = a->row_index[i+1] - a->row_index[i];

        for(usz j = 0; j < N; j++)
        {
            if(row_elems > 0 && a->col_index[idx] == j)
            {
                a_span[i][j] = a->data[idx];
                idx ++;
                row_elems --;
            }
            else
            {
                a_span[i][j] = 0.0;
            }
        }
    }

    return (matrix_t) 
    {
        .size = N,
        .data = new_data
    };
}

matrix_t general_to_general(matrix_t const *a)
{

    const usz N = a->size;
    f64 *new_data = (f64 *)malloc(N * N * sizeof(f64));
    f64 (* new_span)[N] = make_2D_span(f64, , new_data, N);
    f64 (* span)[N] = make_2D_span(f64, , a->data, N);

    for(usz i = 0; i < N; i++)
    {
        for(usz j = 0; j < N; j++)
        {
            new_span[i][j] = span[i][j];
        }
    }

    return (matrix_t) 
    {
        .size = N,
        .data = new_data
    };
}

csr_matrix_t csr_to_csr(csr_matrix_t const *a)
{
    const usz NNZ   = a->row_index[a->size];
    const usz N     = a->size;
    
    f64 *new_data = (f64 *)malloc(NNZ * sizeof(f64));
    usz *new_cols = (usz *)malloc(NNZ * sizeof(usz));
    usz *new_rows = (usz *)malloc((N+1) * sizeof(usz));

    for(usz i = 0; i < NNZ; i++)
    {
        new_data[i] = a->data[i];
        new_cols[i] = a->col_index[i];
    }

    for(usz i = 0; i < N+1; i++)
    {
        new_rows[i] = a->row_index[i];
    }

    return (csr_matrix_t)
    {
        .size = N,
        .row_index  = new_rows,
        .col_index  = new_cols,
        .data       = new_data
    };
}

csr_matrix_t general_to_csr(matrix_t const *a)
{
    const usz NNZ   = NNZ_in_general(a);
    const usz N     = a->size;
    
    f64 *new_data = (f64 *)malloc(NNZ * sizeof(f64));
    usz *new_cols = (usz *)malloc(NNZ * sizeof(usz));
    usz *new_rows = (usz *)malloc((N+1) * sizeof(usz));

    u64 idx = 0;
    new_rows[0] = 0;

    f64 (* a_span)[N] = make_2D_span(f64, , a->data, N);

    for(usz i = 0; i < N; i++)
    {
        for(usz j = 0; j < N; j++)
        { 
            if(a_span[i][j] != 0.0)
            {
                new_data[idx] = a_span[i][j];
                new_cols[idx] = j;
                idx +=1;
            }
        }
        new_rows[i+1] = idx;
    }

    return (csr_matrix_t)
    {
        .size = N,
        .row_index  = new_rows,
        .col_index  = new_cols,
        .data       = new_data
    };
}

u8 csr_equal_csr(csr_matrix_t const *a, csr_matrix_t const *b)
{
    const usz N     = a->size;
    const usz NNZ   = a->row_index[N];
    
    if((N != b->size) && (NNZ != b->row_index[b->size]))
    {
        return 0;
    }
    
    for(usz i = 0; i < NNZ; i++)
    {
        if( (a->data[i] != b->data[i]) && 
            (a->col_index[i] != b->col_index[i]))
        {
            return 0;
        }
    }

    for(usz i = 0; i < N+1; i++)
    {
        if(a->row_index[i] != b->row_index[i])
        {
            return 0;
        }
    }

    return 1;
}

u8 general_equal_csr(matrix_t const *a, csr_matrix_t const *b)
{
    const usz N = a->size;

    if(N != b->size)
    {
        return 0;
    }
    
    u64 idx = 0;
    usz row_elems = 0;
    f64 (* a_span)[N] = make_2D_span(f64, , a->data, N);

    for(usz i = 0; i < N; i++)
    {
        row_elems = b->row_index[i+1] - b->row_index[i];

        for(usz j = 0; j < N; j++)
        {
            if(row_elems > 0 && b->col_index[idx] == j)
            {
                if(a_span[i][j] != b->data[idx])
                {
                    return 0;
                }
                idx ++;
                row_elems --;
            }
            else
            {
               if(a_span[i][j] != 0.0)
               {
                   return 0;
               }
            }
        }
    }

    return 1;
}

u8 general_equal_general(matrix_t const *a, matrix_t const *b)
{ 
    const usz N = a->size;

    if(N != b->size)
    {
        return 0;
    }

    f64 (* a_span)[N] = make_2D_span(f64, , a->data, N);
    f64 (* b_span)[N] = make_2D_span(f64, , b->data, N);

    for(usz i = 0; i < N; i++)
    {
        for(usz j = 0; j < N; j++)
        {
            if(a_span[i][j] != b_span[i][j])
            {
                return 0;
            }
        }
    }

    return 1;
}*/
