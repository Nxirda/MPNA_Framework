#pragma once

#include "lib_matrices.h"
#include "vector.h"

extern void jacobi_general(matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations);
extern void jacobi_csr(csr_matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations);
