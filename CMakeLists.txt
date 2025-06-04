cmake_minimum_required(VERSION 3.21)

if (NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_EXTENSIONS OFF)

project(oledgf
        VERSION 0.1
        LANGUAGES CXX)

include(cmake/PreventInSourceBuilds.cmake)
include(ProjectOptions.cmake)

include(Dependencies.cmake)
oledgf_local_options()
oledgf_setup_dependencies()

find_package(Python REQUIRED Development)
find_package(Python COMPONENTS NumPy)

# Set c++ standard
target_compile_features(oledgf_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_subdirectory(src)
