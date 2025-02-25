#include <stdlib.h>
#include <stdio.h>

#include "solver.h"
#include "fillers.h"
#include "mpi_wrapper.h"

int framework_test(int argc, char **argv)
{
    
    if(argc != 3 && argc != 1)
    {
        fprintf(stderr, "Usage is %s <mesh_size> <nb_iter>\n", argv[0]);
        return 1;
    }
    
    usz mesh_size = 4;
    u64 nb_bytes = mesh_size * mesh_size;

    u64 iter = 1000;
    f64 tol = 1e-6;

    if(argc == 3)
    {
        mesh_size   = strtol(argv[1], NULL, 10);
        nb_bytes    = mesh_size * mesh_size;
        iter        = strtol(argv[2], NULL, 10);
    }

    /* Random target params */
    f64 min = 10e-6;
    f64 max = 10e6;

    vector_t b;
    allocate_vector(&b, nb_bytes);
    
    vector_t x_general;
    allocate_vector(&x_general, nb_bytes);

    vector_t x_csr;
    allocate_vector(&x_csr, nb_bytes);

    csr_matrix_t csr;
    poisson_CSR(mesh_size, &csr);
    
    matrix_t general;
    poisson_general(mesh_size, &general); 
    
    //print_matrix(general);
    /*printf("====== Matrix in general CSC format ======\n");
    csc_matrix_t csc;
    allocate_CSC(mesh_size, &csc);
    fill_CSC(mesh_size, &csc);
    //print_CSC(&csc);*/

    init_random_vector(&b, min, max);
/******************************************************************************/
    {
        printf("====== Jacobi Testing (quick) ======\n");
        
        init_constant_vector(&x_csr, 0.0);
        init_constant_vector(&x_general, 0.0);

        usz gen_iter = jacobi_general(&general, &x_general, &b, iter, tol);
        usz csr_iter = jacobi_csr(&csr, &x_csr, &b, iter, tol);

        u8 jacobi_equal = equal_vector(&x_general, &x_csr); 
        if(jacobi_equal)
            printf("Jacobi implementations yields the same result\n");
        else
            printf("The two implementations yields different results\n");
        
        printf("General iterations : %ld\n", gen_iter);
        printf("Csr iterations     : %ld\n", csr_iter);
        printf("\n");
    }
/******************************************************************************/
    {
        printf("====== Gauss Seidel Testing (quick) ======\n");
        
        init_constant_vector(&x_csr, 0.0);
        init_constant_vector(&x_general, 0.0);

        usz gen_iter = gauss_seidel_general(&general, &x_general, &b, iter, tol);
        usz csr_iter = gauss_seidel_csr(&csr, &x_csr, &b, iter, tol);

        u8 gauss_seidel_equal = equal_vector(&x_general, &x_csr); 
       
        if(gauss_seidel_equal)
            printf("Gauss Seidel implementations yields the same result\n");
        else
            printf("The two implementations yields different results\n");
        
        printf("General iterations : %ld\n", gen_iter);
        printf("Csr iterations     : %ld\n", csr_iter);
        printf("\n");
    }
/******************************************************************************/
    {
        printf("====== Conjugate gradient Testing (quick) ======\n");
        
        init_constant_vector(&x_csr, 0.0);
        init_constant_vector(&x_general, 0.0);
        
        usz gen_iter = conjugate_gradient_general(&general, &x_general, &b, iter, tol);
        usz csr_iter = conjugate_gradient_csr(&csr, &x_csr, &b, iter, tol);

        u8 conjugate_gradient_equal = equal_vector(&x_general, &x_csr); 
        if(conjugate_gradient_equal)
        {
            printf("Conjugate Gradient implementations yields the same result\n");
        }
        else
        {
            printf("The two implementations yields different results\n");
        }
 
        printf("General iterations : %ld\n", gen_iter);
        printf("Csr iterations     : %ld\n", csr_iter);
        printf("\n");
    }
/******************************************************************************/
    {
        printf("====== GMRES Testing (quick) ======\n");
        
        init_constant_vector(&x_csr, 0.0);
        init_constant_vector(&x_general, 0.0);

        usz gen_iter = GMRES_general(&general, &x_general, &b, iter, tol);
        usz csr_iter = GMRES_csr(&csr, &x_csr, &b, iter, tol);
        
        u8 gmres_equal = equal_vector(&x_general, &x_csr); 
        if(gmres_equal)
        {
            printf("Conjugate Gradient implementations yields the same result\n");
        }
        else
        {
            printf("The two implementations yields different results\n");
        }

        printf("General iterations : %ld\n", gen_iter);
        printf("Csr iterations     : %ld\n", csr_iter);
        printf("\n");
    }
    
    free_CSR(&csr);
    free_matrix(general);

    free_vector(&x_general);
    free_vector(&x_csr);
    free_vector(&b);

    return 0;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    i32 size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello from rank %d\n",rank);

    if(argc != 2 && rank == 0)
    {
        printf("Usage is : %s <filename.mtx>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(EXIT_FAILURE);
    }
   
    coo_matrix_t *global_matrix = NULL;
    coo_matrix_t local_matrix;
    vector_t *global_vector = NULL;
    
    if(rank == 0)
    {   
        ascii *filename = argv[1];
        coo_matrix_t mat;
        mm_load_coo(filename, &mat);
        print_COO(&mat);

        global_matrix = &mat;
        
        vector_t vec;
        allocate_vector(&vec, global_matrix->dim_y);
        init_constant_vector(&vec, 2.0);
        global_vector = &vec;
    }
    
    print_COO(global_matrix);
    /*
    distribute_coo(global_matrix, &local_matrix);
    
    csr_matrix_t local_csr;
    coo_to_csr(&local_matrix, &local_csr); 
    
    vector_t local_vector;
    init_random_vector_MPI(&local_vector, local_matrix.dim_x);
    
    print_CSR(&local_csr);
    csr_mv_MPI(&local_csr, &local_vector, global_vector);

    if(rank == 0)
    {
        //print_vector(global_vector);
        free_vector(global_vector);
    }   */
    //print res vec();

    //framework_test(argc, argv);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
