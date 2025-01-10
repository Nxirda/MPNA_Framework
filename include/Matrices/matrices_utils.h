#pragma once

#include "types.h"

#define to_matrix_id(i, j, N) ((i) + N * (j)) 

#define make_2D_span(type, attr, ptr, dim2) (type (attr *)[dim2])(ptr)
