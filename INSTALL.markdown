# IRTK Installation

## Dependencies

1. IRTK requires the following:
     - [FLTK](http://www.fltk.org/) library.
         - Version 1.3 has been tested and is recommended. 
     - [VTK 5 or 6](http://www.vtk.org/) library.
     - [Boost](http://www.boost.org/) library.
         - Not a version earlier than 1.48.
     - The [GNU Scientific Library](http://www.gnu.org/software/gsl/) is also essential.
2. You will also need to have [CMake](http://www.cmake.org/) installed, in order to generate the native build environment.

## Requirements

The applications of the IRTK can be run on Linux distributions
and Mac OS X. On Windows, [Cygwin](https://www.cygwin.com/) must be installed first.

## Build and Installation

1. Clone the repository

    ```shell
    git clone --recurse-submodules https://github.com/BioMedIA/IRTK.git IRTK
    cd IRTK
    ```

2. Create a build directory and change to the build directory

    ```shell
    mkdir build
    cd build
    ```

3. Build IRTK.

    You can configure IRTK with BUILD_TEST option set to ON using [CMake](http://www.cmake.org/cmake/help/runningcmake.html)
    and build the software using the selected build tool (e.g., [GNU Make](http://www.gnu.org/software/make/)) to build the tests.
    ```shell
    cmake ..
    make
    ```

## Configuration

Add the directory containing the IRTK binaries to your **PATH** environment variable.

[GNU Bash](http://www.gnu.org/software/bash/):

```bash
export PATH="$IRTK_DIR/bin:$PATH"
```

[C shell](http://www.computerhope.com/unix/ucsh.htm):

```shell
setenv PATH "$IRTK_DIR/bin:$PATH"
```

where `IRTK_DIR` is the path to the IRTK build or installation directory.
