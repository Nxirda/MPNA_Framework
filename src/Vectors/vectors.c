#include "vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

static inline f64 f64_random(f64 min, f64 max)
{
    assert(min < max);
    f64 res = min;
    f64 dist = max - min;
    f64 frac = (f64)rand()/(f64)RAND_MAX;
    res += dist * frac;
    return res;
}

// We can just flip sign bit tbh which is what __builtin_fabs does 
static inline f64 f64_abs(f64 value)
{
    return value < 0.0 ? -value : value;
}

void allocate_vector(vector_t *vector, usz size)
{
    assert(vector != NULL && size > 0);

    vector->size = size;
    vector->data = (f64 *)malloc(size * sizeof(f64));
}

// For now this is just a sort of std::iota
void init_iota_vector(vector_t *vector, f64 begin)
{
    assert(vector != NULL && vector->data != NULL);

    for(usz i = 0; i < vector->size; i++)
    {
        vector->data[i] = begin + i;
    }
}

void init_random_vector(vector_t *a, f64 min, f64 max)
{ 
    assert(a != NULL && a->data != NULL);

    for(usz i = 0; i < a->size; i++)
    {
        a->data[i] = f64_random(min, max);
    } 
}

void init_constant_vector(vector_t *a, f64 value)
{
    assert(a != NULL && a->data != NULL);

    for(usz i = 0; i < a->size; i++)
    {
        a->data[i] = value;
    }
}

void copy_vector(vector_t const *src, vector_t *dest)
{
    assert(src != NULL && src->data != NULL);
    assert(dest != NULL && dest->data != NULL);
    assert(src->size == dest->size);

    for(usz i = 0; i < src->size; i++)
    {
        dest->data[i] = src->data[i];
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

void swap_vector(vector_t *a, vector_t *b)
{
    assert(a && b && a->data && b->data);

    f64 *tmp = b->data;
    b->data = a->data;
    a->data = tmp;
}


f64 dot_product(vector_t *a, vector_t *b)
{
    const usz N = a->size;
    f64 sum = 0.0;

    for(usz i = 0; i < N; i++)
    {
        sum += a->data[i] * b->data[i];
    }
    return sum;
}

void add_vector(vector_t *a, vector_t *b, vector_t *c)
{
    const usz N = a->size;

    for(usz i = 0; i < N; i++)
    {
        c->data[i] = a->data[i] + b->data[i];   
    }
}

void mul_vector(vector_t *a, vector_t *b, vector_t *c)
{
    const usz N = a->size;

    for(usz i = 0; i < N; i++)
    {
        c->data[i] = a->data[i] * b->data[i];   
    }
}

void add_scalar_vector(f64 scalar, vector_t *a, vector_t *b)
{
    const usz N = a->size;

    for(usz i = 0; i < N; i++)
    {
        b->data[i] = scalar + a->data[i];   
    }
}

void mul_scalar_vector(f64 scalar, vector_t *a, vector_t *b)
{
    const usz N = a->size;

    for(usz i = 0; i < N; i++)
    {
        b->data[i] = scalar * a->data[i];   
    }
}

void daxpy(f64 scalar, vector_t *a, vector_t *b, vector_t *c)
{
    const usz N = a->size;

    for(usz i = 0; i < N; i++)
    {
        c->data[i] = scalar * a->data[i] + b->data[i];   
    }
}

u8 equal_vector(vector_t *a, vector_t *b)
{
    assert(a && b);
    if(a->size != b->size)
    {
        return 0;
    }

    for(usz i = 0; i < a->size; i++)
    {
        if(a->data[i] != b->data[i])
        {
            return 0;
        }
    }

    return 1;
}

// https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
u8 allclose_vector(vector_t *a, vector_t *b, f64 atol, f64 rtol)
{
    assert(a && b);

    if(a->size != b->size)
    {
        return 0;
    }

    f64 lhs = 0.0;
    f64 rhs = 0.0;
    for(usz i = 0; i < a->size; i++)
    { 
        lhs = f64_abs(a->data[i] - b->data[i]);
        rhs = (atol + rtol * f64_abs(b->data[i]));

        if(lhs > rhs)
        {
            return 0;
        }
    }

    return 1;
}
