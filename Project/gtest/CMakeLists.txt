
# gtest/CMakeLists.txt

# Include headers from the main project
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/gtest)

# Add test sources
set(TEST_SOURCES
    runGtest.cpp
   )

# Optionally, you can list header files for clarity
set(TEST_HEADERS
    readFileTest.hpp  # aggiungi qui i file header
)


# Create GTest executable
add_executable(GTest ${TEST_SOURCES} ${TEST_HEADERS})

# Link GTest executable with required libraries
target_link_libraries(GTest PRIVATE ${DFN_LINKED_LIBRARIES})

# Add test
add_test(NAME GTest COMMAND GTest)