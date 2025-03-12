add_library(Source_Lib
    Vectors/vectors.c
    MPI/vector_MPI.c

    Matrices/general.c 
    Matrices/csr.c
    Matrices/csc.c
    
    Matrices/coo.c
    MPI/coo_MPI.c
    Matrices/inter_format.c
    
    Fillers/poisson_laplacian_2D.c
    Fillers/matrix_market_loader.c

    Algorithms/Solvers/jacobi.c
    Algorithms/Solvers/gauss_seidel.c
    Algorithms/Solvers/conjugate_gradient.c
    Algorithms/Solvers/GMRES.c
    Algorithms/Solvers/power_iteration.c

    Algorithms/Linear_Algebra/dgemv.c
    MPI/spmv_MPI.c

    Utils/chrono.c
)

target_link_libraries(Source_Lib
    PUBLIC
    m
    mvec
    LAPACK::LAPACK
    BLAS::BLAS
    lapacke
    MPI::MPI_C
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_include_directories(Source_Lib
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include/Vectors
    ${CMAKE_SOURCE_DIR}/include/Matrices
    ${CMAKE_SOURCE_DIR}/include/Fillers
    ${CMAKE_SOURCE_DIR}/include/MPI
    # Algorithms subfolder, a bit ugly
    ${CMAKE_SOURCE_DIR}/include/Algorithms
    ${CMAKE_SOURCE_DIR}/include/Algorithms/Solvers
    ${CMAKE_SOURCE_DIR}/include/Algorithms/Linear_Algebra
    ${CMAKE_SOURCE_DIR}/include/Utils
)
