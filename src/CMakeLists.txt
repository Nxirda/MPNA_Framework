add_library(Source_Lib
    Vectors/vectors.c

    Matrices/general.c 
    Matrices/csr.c
    Matrices/csc.c
    Matrices/inter_format.c

    Solvers/jacobi.c
    Solvers/gauss_seidel.c

    Utils/chrono.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_include_directories(Source_Lib
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include/Matrices
    ${CMAKE_SOURCE_DIR}/include/Vectors
    ${CMAKE_SOURCE_DIR}/include/Solvers
    ${CMAKE_SOURCE_DIR}/include/Utils
)
