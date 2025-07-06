# oledgf

## Prerequisites

The project is inherently cross-platform and requires:

- CMake
- A C++ compiler
- Python3, Numpy and Matplotlib for plotting

The build system is based on the [cmake_template](https://github.com/cpp-best-practices/cmake_template) by Jason Turner.

For now it includes

- High level warnings and Warnings as Errors (`OFF` by default)
- [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
- [CPM](https://github.com/cpm-cmake/CPM.cmake) for external dependencies

## Getting started

### Clone the GitLab repository

First clone the main repo
```
git clone https://gitlab.ethz.ch/tmarcato/oledgf.git
```
### Build Instructions

By default, the system default compiler will be used. 
CMake uses the environment variable `CXX` to decide which compiler to use.

To configure the build, `ccmake` is encouraged.
From the project root directory you can run
```
ccmake -S . -B ./build
```
Cmake will automatically create the `./build` folder if it does not exist.
Once `ccmake` has finished setting up, press 'c' to configure the project, 'g' to generate the buildsystem and 'q' to quit.

Now you can build the project
```
cd ./build
make 
```

### Contribute
 When contributing and before pushing your changes to a dev branch, please do the following:

 - `clang-format`: run clang-format on all your files using the global .clang-format file from the main branch
