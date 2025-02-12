cmake_minimum_required(VERSION 3.21)

project(oledgf
        VERSION 0.1
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories("${PROJECT_SOURCE_DIR}/include")
add_subdirectory(src)

target_link_libraries(main PUBLIC OLEDFunctions)