#include "test_utils.h"
#include <stdio.h>

void print_test_result(u8 boolean, ascii const *function_name)
{
    if(boolean)
    {
        printf("[O] : Test passed : %s\n", function_name);
    }
    else
    {
        printf("[X] : Test failed : %s\n", function_name);
    }
}
