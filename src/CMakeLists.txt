add_library(OLEDFunctions files.cpp material.cpp solver.cpp)
target_include_directories(OLEDFunctions PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(OLEDFunctions PUBLIC oledgf_options oledgf_warnings PRIVATE Eigen)

add_executable(main main.cpp)

target_link_libraries(main PUBLIC 
                      OLEDFunctions 
                      oledgf_options 
                      oledgf_warnings
                      fmt::fmt
                      Python3::Python
                      Python3::NumPy
                      PRIVATE Eigen)