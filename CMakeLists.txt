cmake_minimum_required (VERSION 2.8)

project (simpleLQCD)
set (simpleLQCD_VERSION_MAJOR 1)
set (simpleLQCD_VERSION_MINOR 0)

set (CMAKE_C_FLAGS "-std=c99")
#set (CMAKE_EXE_LINKER_FLAGS "-lm")

include_directories (./include)
add_subdirectory (modules)
add_executable (main_random ./main/main_random.c)
add_executable (main_ptsource ./main/main_ptsource.c)
target_link_libraries (main_random simplelqcd m)
target_link_libraries (main_ptsource simplelqcd m)
