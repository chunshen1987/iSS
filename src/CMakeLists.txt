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

add_library(${lib_name} SHARED ${SOURCES})
if (APPLE)
    set_target_properties (${lib_name} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (${lib_name} ${LIBS} ${ZLIB_LIBRARIES})

add_executable(${exe_name} main.cpp)
if (APPLE)
    set_target_properties (${exe_name} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (${exe_name} ${lib_name} ${LIBS} ${ZLIB_LIBRARIES})
install(TARGETS ${exe_name} ${lib_name}
        RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)
