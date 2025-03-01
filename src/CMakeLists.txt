add_library(OLEDFunctions files.cpp material.cpp)

add_executable(main main.cpp)

target_include_directories(OLEDFunctions PUBLIC "${PROJECT_SOURCE_DIR}/include")