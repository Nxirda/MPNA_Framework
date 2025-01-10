#include <stdlib.h>
#include <stdio.h>

#include "lib_matrices.h"
#include "jacobi.h"
#include <string.h>

void print_vector_f64(f64 *a, usz size)
{
    for(usz i = 0; i < size; i++)
    {
        printf("%2.3f ", a[i]);
    }
    printf("\n");
}

void iota_vector_f64(f64 *a, f64 begin, usz size)
{
    for(usz i = 0; i < size; i++)
    {
        a[i] = begin;
        begin += 1;
    }
}

int main(int argc, char **argv)
{

    usz mesh_size = 4;
    u64 nb_bytes = mesh_size * mesh_size;
    vector_t b;
    allocate_vector(&b, nb_bytes);
    vector_t x;
    allocate_vector(&x, nb_bytes);

/******************************************************************************/
    csr_matrix_t csr;
    allocate_CSR(mesh_size, &csr);
    fill_CSR(mesh_size, &csr);
    print_CSR(&csr);

    fill_vector(&b, 1.0);
    fill_vector(&x, 3.3);
    
    jacobi_csr(&csr, &x, &b, 1);

    printf("LHS is :\n");
    print_vector(&x);

/******************************************************************************/
    matrix_t general;
    allocate_matrix(mesh_size, &general);
    fill_matrix(mesh_size, general);
    print_matrix(general);
    
    fill_vector(&b, 1.0);
    fill_vector(&x, 3.3);
    
    jacobi_general(&general, &x, &b, 1);

    printf("LHS is :\n");
    print_vector(&x);
/******************************************************************************/ 
    free_CSR(&csr);
    free_matrix(general);
    free_vector(&x);
    free_vector(&b);

    exit(EXIT_SUCCESS);
    return 0;
}
