cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14) #set C++ standard
set(CMAKE_C_STANDARD 11) #set C standard

project(decoding)

#adding dependencies
add_subdirectory(include/mtwister_rng)
add_subdirectory(include/json_library)
add_subdirectory(include/radford_neal)

# Build BP_OSD library
include_directories(src)
set(LIB_FILES ${LIB_FILES} src/bp_decoder_ms.h src/bp_decoder_ms.c)
set(LIB_FILES ${LIB_FILES} src/syndrome.h src/syndrome.c)
set(LIB_FILES ${LIB_FILES} src/binary_char.h src/binary_char.c)
set(LIB_FILES ${LIB_FILES} src/mod2sparse_extra.h src/mod2sparse_extra.c)
set(LIB_FILES ${LIB_FILES} src/load_alist.h src/load_alist.c)
set(LIB_FILES ${LIB_FILES} src/osd.h src/osd.c)
set(LIB_FILES ${LIB_FILES} src/bp_osd.h src/bp_osd.cpp)
set(LIB_FILES ${LIB_FILES} src/sim_util.cpp src/sim_util.h)
set(LIB_FILES ${LIB_FILES} src/timing.cpp src/timing.h)
add_library(bp_osd ${LIB_FILES})
target_link_libraries(bp_osd PUBLIC radford_neal mtwister_rng)
target_link_directories(bp_osd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Build the project files
set(PROJECT bp_osd_decode)
project(${PROJECT})
add_executable(${PROJECT} sim_scripts/bp_osd_decode.cpp)
target_link_libraries(${PROJECT} bp_osd json_library)