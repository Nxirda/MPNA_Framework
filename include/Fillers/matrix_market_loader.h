#pragma once

#include "matrix.h"

extern void mm_load_coo(const ascii *filepath, coo_matrix_t *matrix);
extern void mm_load_csr();
extern void mm_load_csc();
extern void mm_load_general(); 
