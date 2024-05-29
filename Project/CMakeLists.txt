cmake_minimum_required(VERSION 3.9)

project(DFN LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};${CMAKE_CURRENT_SOURCE_DIR}/../externals/Main_Install/eigen3)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};${CMAKE_CURRENT_SOURCE_DIR}/../externals/Main_Install/googletest)

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND DFN_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND DFN_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)

list(APPEND DFN_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND DFN_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND DFN_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND DFN_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources
################################################################################
add_subdirectory(src)
list(APPEND "${CMAKE_NAME_PROJECT}_sources" ${src_sources})
list(APPEND "${CMAKE_NAME_PROJECT}_headers" ${src_headers})
list(APPEND "${CMAKE_NAME_PROJECT}_includes" ${src_includes})
list(APPEND "${CMAKE_NAME_PROJECT}_LINKED_LIBRARIES" ${src_includes})

# Create executable
################################################################################

add_executable(${CMAKE_PROJECT_NAME} main.cpp ${${CMAKE_NAME_PROJECT}_sources} ${${CMAKE_NAME_PROJECT}_headers})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/DFN
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${PROJECT_NAME} ${DFN_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${${CMAKE_NAME_PROJECT}_includes})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create GTest executable
################################################################################
# Assuming you have a `test` directory with test sources like `test_main.cpp`
set(GTEST_SOURCES
	runGtest.cpp  # Aggiungi i tuoi file di test qui
	readFileTest.hpp
)

# Aggiungi tutti gli header file dei test necessari
set(GTEST_HEADERS
    readFileTest.hpp
    # Aggiungi altri file header dei test, se necessario
)

add_executable(runGtest ${GTEST_SOURCES} ${GTEST_HEADERS})
target_link_libraries(runGtest ${DFN_LINKED_LIBRARIES})

# Aggiungi il percorso delle cartelle contenenti i file .hpp
target_include_directories(runGtest PRIVATE 
    ${${CMAKE_NAME_PROJECT}_includes}  # Include le cartelle del progetto principale
    ${CMAKE_SOURCE_DIR}/src  # Include la cartella src del progetto principale
    # Aggiungi altre cartelle di inclusione, se necessario
)

target_compile_options(runGtest PUBLIC -fPIC)