cmake_minimum_required (VERSION 3.4)
project (iSS LANGUAGES CXX)

set(lib_name "iSS")
set(exe_name "iSS.e")

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}") 
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") 

find_package (GSL REQUIRED)
find_package (ZLIB REQUIRED)

if ( ZIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
endif ( ZIB_FOUND)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS} -g -O3 -std=c++11")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -mfpmath=sse")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    find_package(OpenMP)
    if (OpenMP_FOUND)
        set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
    endif ()
endif()

option (link_with_lib "Link executable with the libarary" ON)
option (unittest "Build Unit tests" OFF)

if (unittest)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")
else (unittest)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOCTEST_CONFIG_DISABLE")
endif (unittest)

SET(INCLUDE_DIRS ${GSL_INCLUDE_DIR})
SET(LIBS ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_subdirectory (src)
