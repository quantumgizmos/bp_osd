cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14) #set C++ standard
set(CMAKE_C_STANDARD 11) #set C standard

project(bp_osd_c)

#adding dependencies
add_subdirectory(include/mtwister_rng) #random number generator
add_subdirectory(include/json_library) #JSON library
add_subdirectory(include/radford_neal) #Software for LDPC codes by Radford Neal

# Build BP_OSD library
set(LIB_FILES ${LIB_FILES} src/bp_osd_c.hpp)
set(LIB_FILES ${LIB_FILES} src/bp_decoder_ms.h src/bp_decoder_ms.c)
set(LIB_FILES ${LIB_FILES} src/syndrome.h src/syndrome.c)
set(LIB_FILES ${LIB_FILES} src/binary_char.h src/binary_char.c)
set(LIB_FILES ${LIB_FILES} src/mod2sparse_extra.h src/mod2sparse_extra.c)
set(LIB_FILES ${LIB_FILES} src/load_alist.h src/load_alist.c)
set(LIB_FILES ${LIB_FILES} src/bp_osd.h src/bp_osd.cpp)
set(LIB_FILES ${LIB_FILES} src/sim_util.cpp src/sim_util.h)
set(LIB_FILES ${LIB_FILES} src/timing.cpp src/timing.h)
set(LIB_FILES ${LIB_FILES} src/sort.c src/sort.h)
set(LIB_FILES ${LIB_FILES} src/osd_0.c src/osd_0.h)
set(LIB_FILES ${LIB_FILES} src/osd_w.c src/osd_w.h)
set(LIB_FILES ${LIB_FILES} src/osd_g.c src/osd_g.h)


add_library(bp_osd ${LIB_FILES})
target_link_libraries(bp_osd json_library radford_neal mtwister_rng)
target_include_directories(bp_osd PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src)

## Build the project for the example script
set(PROJECT bp_osd_decode)
project(${PROJECT})
add_executable(${PROJECT} sim_scripts/bp_osd_decode.cpp)
target_link_libraries(${PROJECT} bp_osd)
