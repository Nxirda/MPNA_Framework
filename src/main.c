#include <stdlib.h>
#include <stdio.h>

#include "solver.h"

int main(int argc, char **argv)
{
    
    if(argc != 3 && argc != 1)
    {
        fprintf(stderr, "Usage is %s <mesh_size> <nb_iter>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    
    usz mesh_size = 4;
    u64 nb_bytes = mesh_size * mesh_size;
    u64 iter = 30;
    f64 tol = 1e-6;

    if(argc == 3)
    {
        mesh_size   = strtol(argv[1], NULL, 10);
        nb_bytes    = mesh_size * mesh_size;
        iter        = strtol(argv[2], NULL, 10);
    }

    /* Random target params */
    f64 min = 10e-12;
    f64 max = 10e12;

    vector_t b;
    allocate_vector(&b, nb_bytes);
    
    vector_t x_general;
    allocate_vector(&x_general, nb_bytes);

    vector_t x_csr;
    allocate_vector(&x_csr, nb_bytes);

    printf("====== Matrix in CSR storage format ======\n");
    csr_matrix_t csr;
    allocate_CSR(mesh_size, &csr);
    fill_CSR(mesh_size, &csr);
    //print_CSR(&csr);
   
    printf("====== Matrix in general storage format ======\n");
    matrix_t general;
    allocate_matrix(mesh_size, &general);
    fill_matrix(mesh_size, general);
    //print_matrix(general);
    
    printf("====== Matrix in general CSC format ======\n");
    csc_matrix_t csc;
    allocate_CSC(mesh_size, &csc);
    fill_CSC(mesh_size, &csc);
    //print_CSC(&csc);

    init_random_vector(&b, min, max);
/******************************************************************************/
    printf("====== Jacobi Testing (quick) ======\n");
    
    printf("CSR : \n");

    init_constant_vector(&x_csr, 0.0);
    jacobi_csr(&csr, &x_csr, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_csr);

    printf("General : \n");
    init_constant_vector(&x_general, 0.0);
    jacobi_general(&general, &x_general, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_general);

    u8 jacobi_equal = equal_vector(&x_general, &x_csr); 
    if(jacobi_equal)
    {
        printf("Jacobi implementations yields the same result\n");
    }
    else
    {
        printf("The two implementations yields different results\n");
    }

    printf("\n");
/******************************************************************************/
    printf("====== Gauss Seidel Testing (quick) ======\n");
    
    printf("CSR : \n");
    init_constant_vector(&x_csr, 0.0);
    gauss_seidel_csr(&csr, &x_csr, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_csr);

    printf("General : \n");
    init_constant_vector(&x_general, 0.0);
    gauss_seidel_general(&general, &x_general, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_general);

    u8 gauss_seidel_equal = equal_vector(&x_general, &x_csr); 
    if(jacobi_equal)
    {
        printf("Gauss Seidel implementations yields the same result\n");
    }
    else
    {
        printf("The two implementations yields different results\n");
    }

    printf("\n");
/******************************************************************************/
    printf("====== Conjugate gradient Testing (quick) ======\n");
    
    printf("CSR : \n");
    init_constant_vector(&x_csr, 0.0);
    conjugate_gradient_csr(&csr, &x_csr, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_csr);

    printf("General : \n");
    init_constant_vector(&x_general, 0.0);
    conjugate_gradient_general(&general, &x_general, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_general);

    u8 conjugate_gradient_equal = equal_vector(&x_general, &x_csr); 
    if(conjugate_gradient_equal)
    {
        printf("Conjugate Gradient implementations yields the same result\n");
    }
    else
    {
        printf("The two implementations yields different results\n");
    }

    printf("\n");
/******************************************************************************/
    


    free_CSR(&csr);
    free_matrix(general);

    free_vector(&x_general);
    free_vector(&x_csr);
    free_vector(&b);

    exit(EXIT_SUCCESS);
    //return 0;
}
