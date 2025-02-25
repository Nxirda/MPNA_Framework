cmake_minimum_required(VERSION 3.10)

project(Matrix_Framework LANGUAGES C DESCRIPTION "Framework to experiment with different matrices format and solvers")

set (CMAKE_C_STANDARD 11)
set (CMAKE_C_STANDARD_REQUIRED ON)
set (CMAKE_COLOR_DIAGNOSTICS ON)

set(C_COMPILER "mpicc")

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACKE REQUIRED)

find_package(MPI REQUIRED)

if(C_COMPILER STREQUAL "gcc")
    set(CMAKE_C_COMPILER "gcc")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -g -march=native")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -Wall -Wextra -Wconversions -Wpedantic")
    set(CMAKE_C_LINK_FLAGS_RELEASE "-flto")
    set(CMAKE_C_LINK_FLAGS_DEBUG "")
elseif(C_COMPILER STREQUAL "mpicc")
    set(CMAKE_C_COMPILER "mpicc")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -g -march=native")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -Wall -Wextra -Wconversions -Wpedantic")
    set(CMAKE_C_LINK_FLAGS_RELEASE "-flto")
    set(CMAKE_C_LINK_FLAGS_DEBUG "")
else()
    message(FATAL_ERROR "Invalid C compiler specified: ${C_COMPILER}")
endif()

add_subdirectory(src)
add_executable(main.exe src/main.c)
target_include_directories(main.exe PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(main.exe PRIVATE Source_Lib)

# Tests : Have to ensure compile options are Debug ones here 
enable_testing()

add_subdirectory(test)

add_executable(tests.exe test/main.c)
target_link_libraries(tests.exe PRIVATE Test_lib)
