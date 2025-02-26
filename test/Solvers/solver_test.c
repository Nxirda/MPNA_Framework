#include "solver_test.h"
#include "test_utils.h"

int all_solver_test()
{
    usz mesh_size = 4;
    u64 nb_bytes = mesh_size * mesh_size;

    u64 iter = 1000;
    f64 tol = 1e-6;
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
    
    init_random_vector(&b, min, max);
/******************************************************************************/
    {
        printf("====== Jacobi Testing (quick) ======\n");
        
        init_constant_vector(&x_csr, 0.0);
        init_constant_vector(&x_general, 0.0);

        usz gen_iter = jacobi_general(&general, &x_general, &b, iter, tol);
        usz csr_iter = jacobi_csr(&csr, &x_csr, &b, iter, tol);

        u8 jacobi_equal = equal_vector(&x_general, &x_csr); 
        
        print_test_result(jacobi_equal, "Jacobi Implementation");
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
       
        print_test_result(gauss_seidel_equal, "Gauss Seidel Implementation");
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
        
        print_test_result(conjugate_gradient_equal, "Conjugate Gradient Implementation");
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
        
        print_test_result(gmres_equal, "GEneral Minimal RESidual Implementation");
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

void run_solver_test()
{
    all_solver_test();
}
