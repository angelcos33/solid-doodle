cmake_minimum_required(VERSION 3.21)

if (NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_EXTENSIONS OFF)

project(oledgf
        VERSION 0.1
        LANGUAGES CXX)

include(cmake/PreventInSourceBuilds.cmake)

# Set c++ standard
add_library(oledgf_options INTERFACE)
target_compile_features(oledgf_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_subdirectory(src)

target_link_libraries(main PUBLIC OLEDFunctions oledgf_options)