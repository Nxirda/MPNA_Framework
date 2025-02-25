#include "vector_MPI.h"

void init_random_vector_MPI(vector_t *local_vector, usz vec_size)
{
    i32 size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    allocate_vector(local_vector, vec_size);
    if(rank == 0)
    {
        init_random_vector(local_vector, 10e-7, 10e7);
    }
    //MPI_Bcast(local_vector->data, vec_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
