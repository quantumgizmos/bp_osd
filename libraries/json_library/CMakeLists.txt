add_library(json_library json.hpp json_util.cpp json_util.h)
target_include_directories(json_library PUBLIC ${CMAKE_CURRENT_LIST_DIR})
