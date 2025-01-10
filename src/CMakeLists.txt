add_library(Source_Lib
    Vector/vector.c

    Matrices/general.c    
    Matrices/csr.c

    Solvers/jacobi.c

    Utils/chrono.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_include_directories(Source_Lib
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include/Matrices
    ${CMAKE_SOURCE_DIR}/include/Vector
    ${CMAKE_SOURCE_DIR}/include/Solvers
    ${CMAKE_SOURCE_DIR}/include/Utils
)
