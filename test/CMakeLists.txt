add_Library(
    Test_lib

    Matrices/inter_format_test.c
    Vectors/vector_test.c
    test_utils.c 
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(Test_lib 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/Matrices
    ${CMAKE_CURRENT_SOURCE_DIR}/Vectors
)

target_link_libraries(Test_lib PUBLIC Source_Lib)
