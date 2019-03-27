set (SOURCES
  iSS.cpp
  arsenal.cpp
  ParameterReader.cpp
  RandomVariable1DArray.cpp
  RandomVariable2DArray.cpp
  Random.cpp
  TableFunction.cpp
  Table.cpp
  emissionfunction.cpp
  readindata.cpp
  particle_decay.cpp
  pretty_ostream.cpp
  )

add_library(iSS_lib SHARED ${SOURCES})
if (APPLE)
    set_target_properties (iSS_lib PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (iSS_lib ${LIBS} ${ZLIB_LIBRARIES})

add_executable(iSS.e main.cpp)
if (APPLE)
    set_target_properties (iSS.e PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (iSS.e iSS_lib ${LIBS} ${ZLIB_LIBRARIES})
install(TARGETS iSS.e iSS_lib
        RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)
