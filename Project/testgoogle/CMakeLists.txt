list(APPEND testgoogle_sources "${CMAKE_CURRENT_SOURCE_DIR}/readFileTest.hpp")

list(APPEND testgoogle_headers "${CMAKE_CURRENT_SOURCE_DIR}/readFileTest.hpp")

list(APPEND testgoogle_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(testgoogle_sources ${testgoogle_sources} PARENT_SCOPE)
set(testgoogle_headers ${testgoogle_headers} PARENT_SCOPE)
set(testgoogle_includes ${testgoogle_includes} PARENT_SCOPE)