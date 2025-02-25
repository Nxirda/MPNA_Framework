#pragma once

#include "matrix.h"
#include "vector.h"
#include "algorithms.h"

#include "mpi.h"

extern void csr_mv_MPI(csr_matrix_t *matrix, vector_t *x, vector_t *b);
