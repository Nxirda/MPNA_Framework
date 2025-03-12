#pragma once

#include "matrix.h"
#include "vector.h"

extern usz power_iteration_general(matrix_t const *matrix, vector_t *x, 
                                vector_t const *b, u64 max_iterations);

extern usz power_iteration_csr(csr_matrix_t const *matrix, vector_t *x, 
                                vector_t const *b, u64 max_iterations);
