#include "spmv_MPI.h"
#include <stdio.h>

void csr_mv_MPI(csr_matrix_t *matrix, vector_t *x, vector_t *b)
{
    i32 rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector_t res;
    allocate_vector(&res, x->size);
    init_constant_vector(&res, 0.0);
    
    vector_t fake;
    allocate_vector(&fake, x->size);
    init_constant_vector(&fake, 0.0);
    
    // Mv is actually alpha * A * x + beta * y = res
    csr_mv(1.0, matrix, x, 1.0, &fake, &res);
    //copy_vector(&res, b);
    
    MPI_Allreduce(res.data, b->data, res.size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    free_vector(&res);
    free_vector(&fake);
}
