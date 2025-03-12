#include "power_iteration.h"
#include "linear_algebra.h"
#include <math.h>
#include <float.h> // For floats limits
                   
usz power_iteration_general(matrix_t const *matrix, vector_t *x, 
                            vector_t const *b, u64 max_iterations)
{
    const usz N = x->size;

    vector_t res;
    vector_t bad_api_victim;

    allocate_vector(&res, N);
    allocate_vector(&bad_api_victim, N);
    
    init_constant_vector(&bad_api_victim, 1.0);

    f64 inv_norm = 0.0, lambda = 0.0;

    for(usz it = 0; it < max_iterations; it++)
    {
        inv_norm = dot_product(x, x);
        inv_norm = 1.0/sqrt(inv_norm);

        mul_scalar_vector(inv_norm, x, x);

        general_mv(1.0, matrix, x, 1.0, &bad_api_victim, &res);
        
        swap_vector(x, &res);
    }

    lambda = dot_product(x, &res);

    free_vector(&res);

    return lambda;
}

usz power_iteration_csr(csr_matrix_t const *matrix, vector_t *x, 
                        vector_t const *b, u64 max_iterations)
{
    const usz N = x->size;

    vector_t res;
    vector_t bad_api_victim;

    allocate_vector(&res, N);
    allocate_vector(&bad_api_victim, N);
    
    init_constant_vector(&bad_api_victim, 1.0);

    f64 inv_norm = 0.0, lambda = 0.0;

    for(usz it = 0; it < max_iterations; it++)
    {
        inv_norm = dot_product(x, x);
        inv_norm = 1.0/sqrt(inv_norm);

        mul_scalar_vector(inv_norm, x, x);

        csr_mv(1.0, matrix, x, 1.0, &bad_api_victim, &res);
        
        swap_vector(x, &res);
    }

    lambda = dot_product(x, &res);

    free_vector(&res);

    return lambda;
}
