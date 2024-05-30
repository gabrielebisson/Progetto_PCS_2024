list(APPEND gtest_sources "${CMAKE_CURRENT_SOURCE_DIR}/readFileTest.hpp")

list(APPEND gtest_headers "${CMAKE_CURRENT_SOURCE_DIR}/readFileTest.hpp")

list(APPEND gtest_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gtest_sources ${gtest_sources} PARENT_SCOPE)
set(gtest_headers ${gtest_headers} PARENT_SCOPE)
set(gtest_includes ${gtest_includes} PARENT_SCOPE)