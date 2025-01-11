add_Library(
    Test_lib
    
    inter_format_test.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(Test_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(Test_lib PUBLIC Source_Lib)
