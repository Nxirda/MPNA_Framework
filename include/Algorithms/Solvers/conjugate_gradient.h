#pragma once

#include "matrix.h"
#include "vector.h"

extern usz conjugate_gradient_general(matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol);

extern usz conjugate_gradient_csr(csr_matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations, f64 tol);
