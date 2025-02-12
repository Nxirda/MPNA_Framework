#include "assert.h"

#include "poisson_laplacian_2D.h"

// Specific to poisson
usz laplacian_NNZ_elements(usz mesh_size, usz matrix_size)
{
    assert(mesh_size > 0 && matrix_size > 0);

    const usz other_diagonal_size = matrix_size - mesh_size; 
    return matrix_size + 4 * other_diagonal_size; 
}

void poisson_CSR(usz mesh_size, csr_matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);
    const usz N = mesh_size * mesh_size;
    const usz NNZ = laplacian_NNZ_elements(mesh_size, N);

    allocate_CSR(N, N, NNZ, matrix);

    usz CSR_idx = 0;
    usz row_idx = 0;
    
    matrix->row_index[0] = 0;
    row_idx ++;

    for(usz j = 0; j < mesh_size; j++)
    {
        for(usz i = 0; i < mesh_size; i++) 
        {
            //Far left 
            if(j > 0)
            {
                usz col = to_matrix_id(i,j-1, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }
            //Close left
            if(i > 0)
            {
                usz col = to_matrix_id(i-1,j, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }
           
            usz col = to_matrix_id(i, j, mesh_size);
            matrix->data[CSR_idx] = 4.0;
            matrix->col_index[CSR_idx] = col;
            CSR_idx ++;

            //Close right
            if(i < mesh_size -1)
            {
                usz col = to_matrix_id(i+1,j, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }
            //Far right
            if(j < mesh_size -1)
            {
                usz col = to_matrix_id(i,j+1, mesh_size);
                matrix->data[CSR_idx] = -1.0;
                matrix->col_index[CSR_idx] = col;
                CSR_idx ++;
            }

            matrix->row_index[row_idx] = CSR_idx;
            row_idx ++;
        }
    }
}

void poisson_CSC(usz mesh_size, csc_matrix_t *matrix)
{
    return;
}

void poisson_general(usz mesh_size, matrix_t *matrix)
{
    assert(matrix != NULL && mesh_size > 0);
    const usz N = mesh_size * mesh_size;
    allocate_matrix(N, N, matrix);

    f64 (* span)[matrix->dim_y] = make_2D_span(f64, ,matrix->data, matrix->dim_y);
    
    // This order is due to how we parse the mesh into a matrix
    for(usz j = 0; j < mesh_size; j++)
    {
        for(usz i = 0; i < mesh_size; i++)
        {
            // Corresponding coord in matrix
            const usz k = to_matrix_id(i,j,mesh_size);

            //Far left 
            if(j > 0)
                span[k][to_matrix_id(i,j-1,mesh_size)] = -1.0;
            
            //Close left
            if(i > 0)
                span[k][to_matrix_id(i-1,j,mesh_size)] = -1.0;

            span[k][k] = 4.0;

            //Close right
            if(i < mesh_size -1)
                span[k][to_matrix_id(i+1,j,mesh_size)] = -1.0;
             
            //Far right
            if(j < mesh_size -1)
                span[k][to_matrix_id(i,j+1,mesh_size)] = -1.0;
        }
    }
}
