#pragma once

#include "lib_matrices.h"
#include "vector.h"

void jacobi_general(matrix_t const *matrix, vector_t *x, vector_t const *b, u64 max_iterations);
void jacobi_csr(csr_matrix_t *matrix , f64 *x, f64 *b, usz size);
