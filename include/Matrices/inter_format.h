#pragma once

#include "csr.h"
#include "general.h"

extern usz NNZ_in_general(matrix_t const *a);

extern csr_matrix_t csr_to_csr(csr_matrix_t const *a);
extern csr_matrix_t general_to_csr(matrix_t const *a);
extern matrix_t     csr_to_general(csr_matrix_t const *a);
extern matrix_t     general_to_general(matrix_t const *a);

extern u8 csr_equal_csr(csr_matrix_t const *a, csr_matrix_t const *b);
extern u8 general_equal_csr(matrix_t const *a, csr_matrix_t const *b);
extern u8 general_equal_general(matrix_t const *a, matrix_t const *b);
