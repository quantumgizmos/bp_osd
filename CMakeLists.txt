cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14) #set C++ standard
set(CMAKE_C_STANDARD 11) #set C standard

#direcotry includes
include_directories(mod2sparse_library/radford_neal)
include_directories(bp_osd)
include_directories(random_number_generators/mersenne_twister)
include_directories(json_libary/nlohmann)
include_directories(json_libary)


#Mersenne_twister RNG library load
set(LIB_FILES ${LIB_FILES} random_number_generators/mersenne_twister/mtwister.h random_number_generators/mersenne_twister/mtwister.c)
add_library(mersenne_twister ${LIB_FILES} random_number_generators/mersenne_twister/setup_mtwister.h)

#Build radford neal library
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/mod2sparse.h mod2sparse_library/radford_neal/mod2sparse.c)
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/open.h mod2sparse_library/radford_neal/open.c)
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/mod2dense.h mod2sparse_library/radford_neal/mod2dense.c)
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/mod2convert.h mod2sparse_library/radford_neal/mod2convert.c)
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/intio.h mod2sparse_library/radford_neal/intio.c)
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/rcode.h mod2sparse_library/radford_neal/rcode.c)
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/check.h mod2sparse_library/radford_neal/check.c)
set(LIB_FILES ${LIB_FILES} mod2sparse_library/radford_neal/alloc.h mod2sparse_library/radford_neal/alloc.c)
add_library(radford_neal ${LIB_FILES})

# Build BP_OSD library
set(LIB_FILES ${LIB_FILES} bp_osd/bp_decoder_ms.h bp_osd/bp_decoder_ms.c)
set(LIB_FILES ${LIB_FILES} bp_osd/syndrome.h bp_osd/syndrome.c)
set(LIB_FILES ${LIB_FILES} bp_osd/binary_char.h bp_osd/binary_char.c)
set(LIB_FILES ${LIB_FILES} bp_osd/mod2sparse_extra.h bp_osd/mod2sparse_extra.c)
set(LIB_FILES ${LIB_FILES} bp_osd/load_alist.h bp_osd/load_alist.c)
set(LIB_FILES ${LIB_FILES} bp_osd/osd.h bp_osd/osd.c)
set(LIB_FILES ${LIB_FILES} bp_osd/bp_osd.h bp_osd/bp_osd.cpp)
set(LIB_FILES ${LIB_FILES} bp_osd/sim_util.cpp bp_osd/sim_util.h)
set(LIB_FILES ${LIB_FILES} bp_osd/timing.cpp bp_osd/timing.h)
add_library(bp_osd ${LIB_FILES})

# JSON library
set(LIB_FILES ${LIB_FILES} json_libary/nlohmann/json.hpp json_libary/json_util.h sim_scripts/bp_osd_decode.cpp)
add_library(json_libary ${LIB_FILES})



# Project files
set(PROJECT bp_osd_decode)
project(${PROJECT})
add_executable(${PROJECT} sim_scripts/bp_osd_decode.cpp)
target_link_libraries(${PROJECT} radford_neal)
target_link_libraries(${PROJECT} bp_osd)
target_link_libraries(${PROJECT} mersenne_twister)
target_link_libraries(${PROJECT} json_libary)



# Project files
set(PROJECT joschka_test_boundaries)
project(${PROJECT})
add_executable(${PROJECT} joschka_tests/joschka_test_boundaries.cpp)
target_link_libraries(${PROJECT} radford_neal)
target_link_libraries(${PROJECT} bp_osd)
target_link_libraries(${PROJECT} mersenne_twister)
target_link_libraries(${PROJECT} json_libary)

# Project files
set(PROJECT joschka_test_periodic)
project(${PROJECT})
add_executable(${PROJECT} joschka_tests/joschka_test_periodic.cpp)
target_link_libraries(${PROJECT} radford_neal)
target_link_libraries(${PROJECT} bp_osd)
target_link_libraries(${PROJECT} mersenne_twister)
target_link_libraries(${PROJECT} json_libary)


# Project files
set(PROJECT joschka_test_periodic2)
project(${PROJECT})
add_executable(${PROJECT} joschka_tests/joschka_test_periodic2.cpp)
target_link_libraries(${PROJECT} radford_neal)
target_link_libraries(${PROJECT} bp_osd)
target_link_libraries(${PROJECT} mersenne_twister)
target_link_libraries(${PROJECT} json_libary)

# Project files
set(PROJECT joschka_test_periodic3)
project(${PROJECT})
add_executable(${PROJECT} joschka_tests/joschka_test_periodic3.cpp)
target_link_libraries(${PROJECT} radford_neal)
target_link_libraries(${PROJECT} bp_osd)
target_link_libraries(${PROJECT} mersenne_twister)
target_link_libraries(${PROJECT} json_libary)