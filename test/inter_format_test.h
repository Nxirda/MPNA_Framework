#pragma once

#include "matrix.h"

/* Helper function */
extern void run_inter_format_test();
static void print_test_result(u8 boolean, ascii const *function_name);

/**/
static void test_NNZ_in_general_count();

/* Equality functions between two formats */
static void test_general_equal_general();
static void test_general_equal_csr();
static void test_csr_equal_csr();

/* Copy functions from base format to target format */
static void test_csr_to_general();
static void test_csr_to_csr();
static void test_general_to_general();
static void test_general_to_csr();
