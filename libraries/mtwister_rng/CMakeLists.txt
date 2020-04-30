add_library(mtwister_rng mtwister.h mtwister.c setup_mtwister.h setup_mtwister.cpp)
target_include_directories(mtwister_rng PUBLIC ${CMAKE_CURRENT_LIST_DIR})
