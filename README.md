# oledgf

## Prerequisites

### GCC & Make

- Linux:
    - Debian/Ubuntu: `sudo apt install build-essential`
    - Fedora: `sudo dnf install gcc-c++ make`
    - Arch: `sudo pacman -S base-devel`
- macOS:
    1. Run the following command: `xcode-select --install`
    2. In the window which pops up, click "Install" and follow the instructions.
- Windows:
    1. Install MinGW-w64 from [WinLibs.com](https://winlibs.com/).
    2. Add the path to MinGW-64's `bin` directory to Windows's system `PATH` environment variable.
        > You will need to use `mingw32-make` instead of `make` each time the `make` command is used in this README.
    3. Install Git Bash by installing [Git for Windows](https://git-scm.com/downloads).
        > You will need to use **Git Bash** over PowerShell or cmd.exe each time the `make` command is used in this README.


## Build Application

To build the application you first have to clone the repo

```
git clone https://gitlab.ethz.ch/tmarcato/oledgf.git
```

 Then you can build it using the provided Makefile
```
make 
```

### Contribute
 When contributing and before pushing your changes to a dev branch, please do the following:

 - `clang-format`: run clang-format on all your files using the global .clang-format file from the main branch
 - `make clean`
