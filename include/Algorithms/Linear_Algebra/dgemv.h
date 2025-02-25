#pragma once

#include "matrix.h"
#include "vector.h"

extern void general_mv(f64 alpha, matrix_t *A, vector_t *x, f64 beta, vector_t *b, vector_t *r);
extern void csr_mv(f64 alpha, csr_matrix_t *A, vector_t *x, f64 beta, vector_t *b, vector_t *r);
