# how to write CMakeList.txt: https://subingwen.cn/cmake/CMake-primer
cmake_minimum_required (VERSION 3.12.0)

# Project information
project (BaseVar VERSION 1.1.0 LANGUAGES CXX)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# find all source codes
file(GLOB SRC_LIST ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/htslib)      # include the header of htslib
link_libraries(${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)  # link htslib static lib: htslib.a

# create executable basevar
add_executable(basevar ${SRC_LIST})

# Compile htslib
add_custom_target(
    Buildhts ALL
    COMMAND autoreconf -i && ./configure && make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
)

set(COMMON_SYS_LIB pthread z bz2 m lzma curl)
if (APPLE)
    target_link_libraries(basevar ${COMMON_SYS_LIB})
else ()
    target_link_libraries(basevar ${COMMON_SYS_LIB} ssl crypto)
endif()
 
