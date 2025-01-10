#include "vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void allocate_vector(vector_t *vector, usz size)
{
    assert(vector != NULL && size > 0);

    vector->size = size;
    vector->data = (f64 *)malloc(size * sizeof(f64));
}

// For now this is just a sort of std::iota
void fill_vector(vector_t *vector, usz begin)
{
    assert(vector != NULL && vector->data != NULL);

    for(usz i = 0; i < vector->size; i++)
    {
        vector->data[i] = begin + i;
    }
}

void print_vector(vector_t *vector)
{
    assert(vector != NULL && vector->data != NULL);

    for(usz i = 0; i < vector->size; i++)
    {
        printf("%2.3f ", vector->data[i]);
    }
    printf("\n");
}

void free_vector(vector_t *vector)
{
    assert(vector != NULL && vector->data != NULL);
    free(vector->data);
}
