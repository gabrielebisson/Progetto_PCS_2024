list(APPEND visualPara_sources "${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.hpp")
list(APPEND visualPara_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils2.hpp")
list(APPEND visualPara_sources "${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.cpp")
list(APPEND visualPara_sources "${CMAKE_CURRENT_SOURCE_DIR}/ExportDFN.cpp")

list(APPEND visualPara_headers "${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.hpp")
list(APPEND visualPara_headers "${CMAKE_CURRENT_SOURCE_DIR}/Utils2.hpp")

list(APPEND visualPara_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(visualPara_sources ${visualPara_sources} PARENT_SCOPE)
set(visualPara_headers ${visualPara_headers} PARENT_SCOPE)
set(visualPara_includes ${visualPara_includes} PARENT_SCOPE)