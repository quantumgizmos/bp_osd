# single-shot
Single-shot decoding of topological codes using MWPM &amp; BP+OSD

## Build instructions

Use [CMake](https://cmake.org/) to build. During the build process CMake will automatically download [googletest](https://github.com/google/googletest).

### Step-by-step instructions

- `mkdir build && cd build`
- `cmake -DCMAKE_BUILD_TYPE=Release ../`
- `make`

### To run the tests

- `mkdir build && cd build`
- `cmake -DCMAKE_BUILD_TYPE=Debug ../`
- `make`
- `make test`

Tested on Linux (Ubuntu 16.04).

## BP OSD Instructions

### Mod2libary setup
Currently this has to be compiled manually. There must be a way of doing this automatically from CMAKE but I couldn't figure it out...

- `cd mod2sparse_library/radford_neal`
- `make`

### To run bp_osd_hamming.cpp

- Change working directory to the repository root.
- Run "joschka" project.
