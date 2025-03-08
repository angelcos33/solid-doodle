add_library(OLEDFunctions files.cpp material.cpp)
target_include_directories(OLEDFunctions PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(OLEDFunctions PUBLIC oledgf_options)

add_executable(main main.cpp)
