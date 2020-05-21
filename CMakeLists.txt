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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -mfpmath=sse")

option (link_with_lib "Link executable with the libarary" ON)

SET(INCLUDE_DIRS ${GSL_INCLUDE_DIR})
SET(LIBS ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_subdirectory (src)
