#include <stdlib.h>
#include <stdio.h>

#include "solver.h"

int main(int argc, char **argv)
{

    usz mesh_size = 4;
    u64 nb_bytes = mesh_size * mesh_size;
    u64 iter = 30;

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
    print_CSR(&csr);
   
    printf("====== Matrix in general storage format ======\n");
    matrix_t general;
    allocate_matrix(mesh_size, &general);
    fill_matrix(mesh_size, general);
    print_matrix(general);

    init_random_vector(&b, min, max);
/******************************************************************************/
    printf("====== Jacobi Testing (quick) ======\n");
    
    printf("CSR : \n");

    init_constant_vector(&x_csr, 0.0);
    jacobi_csr(&csr, &x_csr, &b, iter);
    printf("LHS is :\n");
    print_vector(&x_csr);

    printf("General : \n");
    init_constant_vector(&x_general, 0.0);
    jacobi_general(&general, &x_general, &b, iter);
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
    gauss_seidel_csr(&csr, &x_csr, &b, iter);
    printf("LHS is :\n");
    print_vector(&x_csr);

    printf("General : \n");
    init_constant_vector(&x_general, 0.0);
    gauss_seidel_general(&general, &x_general, &b, iter);
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

    free_CSR(&csr);
    free_matrix(general);

    free_vector(&x_general);
    free_vector(&x_csr);
    free_vector(&b);

    exit(EXIT_SUCCESS);
    return 0;
}
