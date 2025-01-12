#include "vector_test.h"
#include "test_utils.h"

const usz size = 15;

void test_init_iota()
{
    f64 begin = 0.0;
    u8 ret = 1;
    vector_t a;
    allocate_vector(&a, size);

    init_iota_vector(&a, begin);
    
    for(usz i = 0; i < size; i++)
    {
        if(a.data[i] != begin + i)
        {
            ret = 0;
        }
    }

    print_test_result(ret, "vector iota init");
    free_vector(&a);
}

void test_init_random()
{
    f64 min = 10e-13;
    f64 max = 10e13; 
    u8 ret = 1;
    vector_t a;
    allocate_vector(&a, size);

    init_random_vector(&a, min, max);
    
    for(usz i = 0; i < size; i++)
    {
        if((a.data[i] > max) || (a.data[i] < min))
        {
            ret = 0;
        }
    }

    print_test_result(ret, "vector random init");
    free_vector(&a);
}

void test_init_constant()
{
    f64 value = 456123.312638123; 
    u8 ret = 1;
    vector_t a;
    allocate_vector(&a, size);

    init_constant_vector(&a, value);
    
    for(usz i = 0; i < size; i++)
    {
        if(a.data[i] != value)
        {
            ret = 0;
        }
    }

    print_test_result(ret, "vector constant init");
    free_vector(&a);
}

void test_swap_vector()
{
    f64 value = 456123.312638123; 
    u8 ret = 1;

    vector_t a;
    vector_t b;
    
    allocate_vector(&a, size);
    allocate_vector(&b, size);
    
    init_constant_vector(&a, value);
    init_constant_vector(&b, value/2.0);

    swap_vector(&a, &b);

    for(usz i = 0; i < size; i++)
    {
        if((b.data[i] != value) || (a.data[i] != value/2.0))
        {
            ret = 0;
        }
    }

    print_test_result(ret, "vector swapping");
    free_vector(&a);
    free_vector(&b);
}

void test_equal_vector()
{
    f64 value = 42.69;
    u8 ret;
    vector_t a;
    vector_t b;
  
    allocate_vector(&a, size);
    allocate_vector(&b, size);

    init_constant_vector(&a, value);
    init_constant_vector(&b, value);

    ret = equal_vector(&a, &b);
    print_test_result(ret, "vector equal (same vectors)");

    init_constant_vector(&a, value);
    init_constant_vector(&b, value +1);

    ret = equal_vector(&a, &b);
    print_test_result(!ret, "vector equal (different vectors)");

    free_vector(&a);
    free_vector(&b);
}

void test_allclose_vector()
{
    f64 value = 42.69;
    f64 atol = 10e-8;
    f64 rtol = 10e-5;

    u8 ret;
    vector_t a;
    vector_t b;
    
    init_constant_vector(&a, value);
    init_constant_vector(&b, value);

    ret = allclose_vector(&a, &b, atol, rtol);
    print_test_result(ret, "vector all close (same vectors)");

    init_constant_vector(&a, value);
    init_constant_vector(&b, value + 10e-13);

    ret = allclose_vector(&a, &b, atol, rtol);
    print_test_result(ret, "vector all close (different vectors : variation outside threshold)");

    init_constant_vector(&a, value);
    init_constant_vector(&b, value + 10e-2);

    ret = allclose_vector(&a, &b, atol, rtol);
    print_test_result(!ret, "vector all close (different vectors : variation inside threshold)");

    free_vector(&a);
    free_vector(&b);
}

void run_vector_test()
{
    test_init_iota();
    test_init_random();
    test_init_constant();
    test_equal_vector();
    test_allclose_vector();
}
