cmake_minimum_required (VERSION 2.6)
project (iSS)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

find_package (GSL REQUIRED)
find_package (ZLIB REQUIRED)

if ( ZIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
endif ( ZIB_FOUND)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS} -g -O3 -Wall")

SET(INCLUDE_DIRS ${GSL_INCLUDE_DIR})
SET(LIBS ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_subdirectory (src)
