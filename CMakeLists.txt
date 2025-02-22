cmake_minimum_required(VERSION 3.21)

if (NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

project(oledgf
        VERSION 0.1
        LANGUAGES CXX)

include(cmake/PreventInSourceBuilds.cmake)

include_directories("${PROJECT_SOURCE_DIR}/include")
add_subdirectory(src)

target_link_libraries(main PUBLIC OLEDFunctions)