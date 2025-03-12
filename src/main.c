#include <stdlib.h>
#include <stdio.h>

#include "solver.h"
#include "fillers.h"
#include "mpi_wrapper.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    i32 size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc != 2 && rank == 0)
    {
        printf("Usage is : %s <filename.mtx>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(EXIT_FAILURE);
    }
   
    coo_matrix_t *global_matrix = NULL;
    coo_matrix_t local_matrix;
    
    if(rank == 0)
    {   
        ascii *filename = argv[1];
        coo_matrix_t mat;
        mm_load_coo(filename, &mat);
        global_matrix = &mat;
    }
    
    distribute_coo(global_matrix, &local_matrix);

    csr_matrix_t local_csr;
    coo_to_csr(&local_matrix, &local_csr); 
    
    vector_t local_vector;
    init_random_vector_MPI(&local_vector, local_matrix.dim_x);
   
    vector_t result;
    allocate_vector(&result, local_csr.dim_y);
    
    csr_mv_MPI(&local_csr, &local_vector, &result);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
    {
        csr_matrix_t sequential_mat;
        coo_to_csr(global_matrix, &sequential_mat);

        vector_t sequential_res;
        allocate_vector(&sequential_res, sequential_mat.dim_y);
        
        vector_t fake_b;
        allocate_vector(&fake_b, local_vector.size);
        init_constant_vector(&fake_b, 0.0);

        // Actually performs a alpha * A * x + beta * b = res operation ...
        // I just tweaked the code and didnt take time to rename yet
        csr_mv(1.0, &sequential_mat, &local_vector, 1.0, &fake_b, &sequential_res);
        
        u8 res = equal_vector(&sequential_res, &result);
        
        if(res)
        {
            printf("Sequential and Distributed sparse matrix vector implementations yields the same result\n");
        }
        else
        {
            printf("The two implementations of sparse matrix vector yields different results\n");
        }
        free_vector(&sequential_res);
        free_CSR(&sequential_mat);
        free(&global_matrix);
    }

    free(&local_vector);
    free_CSR(&local_csr);
    free_vector(&result);

    free_COO(&local_matrix);
    //framework_test(argc, argv);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
