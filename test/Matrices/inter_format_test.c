#include "inter_format.h"

#include "test_utils.h"
#include <stdio.h>

const usz mesh_size = 4;

void test_NNZ_in_general_count()
{
    matrix_t a;
    csr_matrix_t b;

    allocate_matrix(mesh_size, &a);
    allocate_CSR(mesh_size, &b);

    // Filling two poisson matrices of the same size
    // so NNZ elems shall be the same 
    fill_matrix(mesh_size, a);
    fill_CSR(mesh_size, &b);

    usz NNZ_a = NNZ_in_general(&a);
    usz NNZ_b = b.row_index[b.size];

    u8 boolean = NNZ_a == NNZ_b;
    
    print_test_result(boolean, "NNZ count in general");
    
    free_matrix(a);
    free_CSR(&b);
}

void test_general_equal_general()
{
    matrix_t a;
    matrix_t b;
    
    allocate_matrix(mesh_size, &a);
    allocate_matrix(mesh_size, &b);

    fill_matrix(mesh_size, a);
    fill_matrix(mesh_size, b);

    u8 geg = general_equal_general(&a, &b);
    print_test_result(geg, "general to general equality");

    free_matrix(a);
    free_matrix(b);
}

void test_general_equal_csr()
{
    matrix_t a;
    csr_matrix_t b;
    
    allocate_matrix(mesh_size, &a);
    allocate_CSR(mesh_size, &b);

    fill_matrix(mesh_size, a);
    fill_CSR(mesh_size, &b);

    u8 gec = general_equal_csr(&a, &b);
    print_test_result(gec, "general to csr equality");

    free_matrix(a);
    free_CSR(&b);
}

void test_csr_equal_csr()
{
    csr_matrix_t a;
    csr_matrix_t b;
    
    allocate_CSR(mesh_size, &a);
    allocate_CSR(mesh_size, &b);

    fill_CSR(mesh_size, &a);
    fill_CSR(mesh_size, &b);
    
    u8 cec = csr_equal_csr(&a, &b);
    print_test_result(cec, "csr to csr equality");

    free_CSR(&a);
    free_CSR(&b);
}

void test_csr_to_general()
{
    csr_matrix_t base;
    allocate_CSR(mesh_size, &base);
    fill_CSR(mesh_size, &base);
    matrix_t copy = csr_to_general(&base);
    
    u8 gtc = general_equal_csr(&copy, &base);
    print_test_result(gtc, "csr to general matrix copy");

    free_CSR(&base);
    free_matrix(copy);
}

void test_csr_to_csr()
{
    csr_matrix_t base;
    allocate_CSR(mesh_size, &base);
    fill_CSR(mesh_size, &base);
    
    csr_matrix_t copy = csr_to_csr(&base);
    u8 ctc = csr_equal_csr(&copy, &base); 
    print_test_result(ctc, "csr to csr matrix copy");

    free_CSR(&base);
    free_CSR(&copy);
}

void test_general_to_general()
{
    matrix_t base;
    allocate_matrix(mesh_size, &base);
    fill_matrix(mesh_size, base);
    
    matrix_t copy = general_to_general(&base);
    u8 gtg = general_equal_general(&copy, &base);
    print_test_result(gtg, "general to general matrix copy");

    free_matrix(base);
    free_matrix(copy);
}

void test_general_to_csr()
{
    matrix_t base;
    allocate_matrix(mesh_size, &base);
    fill_matrix(mesh_size, base);
    
    csr_matrix_t copy = general_to_csr(&base);

    u8 gtc = general_equal_csr(&base, &copy);
    print_test_result(gtc, "general to csr matrix copy");

    free_matrix(base);
    free_CSR(&copy);
}

void run_inter_format_test()
{

    test_NNZ_in_general_count();

    test_general_equal_general();
    test_general_equal_csr();
    test_csr_equal_csr();
    
    test_csr_to_general();
    test_csr_to_csr();
    test_general_to_general();
    test_general_to_csr();
}
