#include <stdlib.h>
#include <stdio.h>

#include "solver.h"
#include "poisson_laplacian_2D.h"

int framework_test(int argc, char **argv)
{
    
    if(argc != 3 && argc != 1)
    {
        fprintf(stderr, "Usage is %s <mesh_size> <nb_iter>\n", argv[0]);
        exit(EXIT_SUCCESS);
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

    //printf("====== Matrix in CSR storage format ======\n");
    csr_matrix_t csr;
    poisson_CSR(mesh_size, &csr);
    //allocate_CSR(mesh_size, &csr);
    //fill_CSR(mesh_size, &csr);
    //print_CSR(&csr);
   
    //printf("====== Matrix in general storage format ======\n");
    matrix_t general;
    poisson_general(mesh_size, &general); 
    //allocate_matrix(mesh_size, &general);
    //fill_matrix(mesh_size, general);
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
        usz jc_iter = jacobi_csr(&csr, &x_csr, &b, iter, tol);

        init_constant_vector(&x_general, 0.0);
        usz jg_iter = jacobi_general(&general, &x_general, &b, iter, tol);

        u8 jacobi_equal = equal_vector(&x_general, &x_csr); 
        if(jacobi_equal)
            printf("Jacobi implementations yields the same result\n");
        else
            printf("The two implementations yields different results\n");
        
        printf("General iterations : %ld\n", jg_iter);
        printf("Csr iterations     : %ld\n", jc_iter);
        printf("\n");
    }
/******************************************************************************/
    {
        printf("====== Gauss Seidel Testing (quick) ======\n");
        
        init_constant_vector(&x_csr, 0.0);
        init_constant_vector(&x_general, 0.0);

        usz gsc_iter = gauss_seidel_csr(&csr, &x_csr, &b, iter, tol);
        usz gsg_iter = gauss_seidel_general(&general, &x_general, &b, iter, tol);

        u8 gauss_seidel_equal = equal_vector(&x_general, &x_csr); 
       
        if(gauss_seidel_equal)
            printf("Gauss Seidel implementations yields the same result\n");
        else
            printf("The two implementations yields different results\n");
        
        printf("General iterations : %ld\n", gsg_iter);
        printf("Csr iterations     : %ld\n", gsc_iter);
        printf("\n");
    }
/******************************************************************************/
    {
        printf("====== Conjugate gradient Testing (quick) ======\n");
        
        //printf("CSR : \n");
        init_constant_vector(&x_csr, 0.0);
        usz cjc_iter = conjugate_gradient_csr(&csr, &x_csr, &b, iter, tol);
        //printf("LHS is :\n");
        //print_vector(&x_csr);

        //printf("General : \n");
        init_constant_vector(&x_general, 0.0);
        usz cjg_iter = conjugate_gradient_general(&general, &x_general, &b, iter, tol);
        //printf("LHS is :\n");
        //print_vector(&x_general);

        u8 conjugate_gradient_equal = equal_vector(&x_general, &x_csr); 
        if(conjugate_gradient_equal)
        {
            printf("Conjugate Gradient implementations yields the same result\n");
        }
        else
        {
            printf("The two implementations yields different results\n");
        }
 
        printf("General iterations : %ld\n", cjg_iter);
        printf("Csr iterations     : %ld\n", cjc_iter);
        printf("\n");
    }
/******************************************************************************/
    /*printf("====== GMRES Testing (quick) ======\n");
    
    printf("CSR : \n");
    init_constant_vector(&x_csr, 0.0);
    conjugate_gradient_csr(&csr, &x_csr, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_csr);

    printf("General : \n");
    init_constant_vector(&x_general, 0.0);
    GMRES_general(&general, &x_general, &b, iter, tol);
    printf("LHS is :\n");
    print_vector(&x_general);
*/
    /*u8 conjugate_gradient_equal = equal_vector(&x_general, &x_csr); 
    if(conjugate_gradient_equal)
    {
        printf("Conjugate Gradient implementations yields the same result\n");
    }
    else
    {
        printf("The two implementations yields different results\n");
    }*/

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

#define MAX_SIZE 2048
void init_Matrix_Market_CSR(ascii *filename, csr_matrix_t *matrix)
{
    FILE *file = fopen(filename, "r");
    if(!file)
    {
        perror("Failed to open given file");
        exit(EXIT_FAILURE);
    }
    
    ascii buffer[MAX_SIZE];
    while(fgets(buffer, sizeof(buffer), file)) 
    {
        if(buffer[0] == '%') continue;
        break;
    }

    usz rows, cols, NNZ;
    sscanf(buffer, "%ld %ld %ld", &rows, &cols, &NNZ);

    matrix->data        = (f64 *)malloc(NNZ * sizeof(f64));
    matrix->col_index   = (usz *)malloc(NNZ * sizeof(usz));
    matrix->row_index   = (usz *)malloc((rows+1) * sizeof(usz));
    matrix->size        = rows * cols;
    
    usz row = 0;
    usz col = 0;
    usz curr_nnz = 0;
    f64 value = 0.0;
    while(fgets(buffer, sizeof(buffer), file))
    {
        if(buffer[0] == '%') continue;

        sscanf(buffer, "%ld %ld %lf", &row, &col, &value);
      
        // As the format is a 1 based indexing
        row --;
        col --;

        matrix->data[curr_nnz] = value;
        matrix->col_index[curr_nnz] = col;
        curr_nnz ++;
    
        matrix->row_index[row+1]++;
    }
    // Convert row_ptr to cumulative sum
    for (int i = 1; i <= rows; i++) {
        matrix->row_index[i] += matrix->row_index[i - 1];
    }
    fclose(file);
}

//void mpi_gemv_csr(csr_matrix_t const* base, usz nb_procs)

int main(int argc, char **argv)
{
    /*if(argc != 2)
    {
        printf("Usage is : %s <filename.mtx>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    // Read in COO maybe ?
    // Read and fill matrix on rank 0
    // Distribute on all ranks

    ascii *filename = argv[1];

    csr_matrix_t mat;
    init_Matrix_Market_CSR(filename, &mat);
    print_CSR(&mat);*/
    
    /*vector_t vec;
    allocate_vector(&vec, mat.size);
    init_constant_vector(&vec, 0.0);

    vector_t res;
    allocate_vector(&res, mat.size);
    init_constant_vector(&res, 0.0);*/
    
    //gemv_CSR()
    
    framework_test(argc, argv);
    exit(EXIT_SUCCESS);
}
