#pragma once

#include "matrix.h"
#include <mpi.h>

extern void distribute_coo(coo_matrix_t *global_matrix, coo_matrix_t *local_matrix);
