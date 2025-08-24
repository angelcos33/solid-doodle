cmake_minimum_required(VERSION 3.29)

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

set (CMAKE_FIND_FRAMEWORK NEVER)
set (Python3_ARTIFACTS_INTERACTIVE TRUE)
set (Python3_FIND_VIRTUALENV FIRST)
set (Python3_FIND_STRATEGY LOCATION)
find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter Development.Embed)
message(STATUS "${Python3_FOUND}")
message(STATUS "${Python3_EXECUTABLE}")
message(STATUS "${Python3_SITEARCH}")
message(STATUS "${Python3_Development.Embed_FOUND}")
message(STATUS "${Python3_INCLUDE_DIRS}")
message(STATUS "${Python3_LIBRARIES}")
find_package(Python3 COMPONENTS NumPy)
message(STATUS "${Python3_Numpy_FOUND}")
message(STATUS "${Python3_Numpy_INCLUDE_DIRS}")
set(CMAKE_FIND_FRAMEWORK FIRST)

# Set c++ standard
target_compile_features(oledgf_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_subdirectory(src)
