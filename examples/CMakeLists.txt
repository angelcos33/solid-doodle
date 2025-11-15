add_executable(test_fitting test_fitting.cpp)
target_link_libraries(test_fitting
                      PUBLIC oledgf_warnings
                             oledgf_options
                             matplot
                      PRIVATE Eigen)               