cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14) #set C++ standard
set(CMAKE_C_STANDARD 11) #set C standard

#direcotry includes
include_directories(src)
include_directories(libraries/mod2sparse_library/radford_neal)
include_directories(libraries/random_number_generators/mersenne_twister)
include_directories(libraries/json_libary/nlohmann)
include_directories(libraries/json_libary)


#Mersenne_twister RNG library load
set(LIB_FILES ${LIB_FILES} libraries/random_number_generators/mersenne_twister/mtwister.h libraries/random_number_generators/mersenne_twister/mtwister.c)
add_library(mersenne_twister ${LIB_FILES} libraries/random_number_generators/mersenne_twister/setup_mtwister.h)

#Build radford neal library
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/mod2sparse.h libraries/mod2sparse_library/radford_neal/mod2sparse.c)
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/open.h libraries/mod2sparse_library/radford_neal/open.c)
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/mod2dense.h libraries/mod2sparse_library/radford_neal/mod2dense.c)
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/mod2convert.h libraries/mod2sparse_library/radford_neal/mod2convert.c)
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/intio.h libraries/mod2sparse_library/radford_neal/intio.c)
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/rcode.h libraries/mod2sparse_library/radford_neal/rcode.c)
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/check.h libraries/mod2sparse_library/radford_neal/check.c)
set(LIB_FILES ${LIB_FILES} libraries/mod2sparse_library/radford_neal/alloc.h libraries/mod2sparse_library/radford_neal/alloc.c)
add_library(radford_neal ${LIB_FILES})

# Build BP_OSD library
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

# JSON library
set(LIB_FILES ${LIB_FILES} libraries/json_libary/nlohmann/json.hpp libraries/json_libary/json_util.h sim_scripts/bp_osd_decode.cpp)
add_library(json_libary ${LIB_FILES})



# Project files
set(PROJECT bp_osd_decode)
project(${PROJECT})
add_executable(${PROJECT} sim_scripts/bp_osd_decode.cpp)
target_link_libraries(${PROJECT} radford_neal)
target_link_libraries(${PROJECT} bp_osd)
target_link_libraries(${PROJECT} mersenne_twister)
target_link_libraries(${PROJECT} json_libary)

