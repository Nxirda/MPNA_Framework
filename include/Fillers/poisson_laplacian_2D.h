#pragma once

#include "matrix.h"

extern void poisson_CSR(usz mesh_size, csr_matrix_t *matrix);
extern void poisson_CSC(usz mesh_size, csc_matrix_t *matrix);
extern void poisson_general(usz mesh_size, matrix_t *matrix);
