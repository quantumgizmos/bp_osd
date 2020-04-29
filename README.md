# BP+OSD (C version)
A C library implementing belief propagation with ordered statistics post-processing for decoding quantum codes.

## Build
To build, run the following commands from the repository root.

```
mkdir build
cd build
cmake ../
make
cd ..
```

## Run BP+OSD Sim
The C++ source file for the BP+OSD sim script can be found at sim_scripts/bp_osd_decode.cpp
The example program to decode a 16_4_6 MKMN code can be run (from the repository root) with the following command

```
build/bp_osd_decode sim_scripts/mkmn_input_test.json sim_scripts/output
```

The general syntax for running the sim script is as follows

```
./build/bp_osd_decode <path_to_input_file.json> <path_to_output_directory>
``` 
