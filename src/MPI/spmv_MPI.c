#include "spmv_MPI.h"
#include <stdio.h>

void csr_mv_MPI(csr_matrix_t *matrix, vector_t *x, vector_t *b)
{
    printf("here\n");
    i32 rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector_t res;
    allocate_vector(&res, x->size);
    init_constant_vector(&res, 0.0);
    
    vector_t fake;
    allocate_vector(&fake, x->size);
    init_constant_vector(&fake, 0.0);
    
    print_vector(&res);
    printf("\n");
    print_vector(&fake);

    csr_mv(1.0, matrix, x, 1.0, &fake, &res);
    copy_vector(&res, b);
    print_vector(&res);
    //MPI_Allreduce(res.data, b->data, x->size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
