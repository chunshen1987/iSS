set (SOURCES
  main.cpp 
  iSS.cpp
  arsenal.cpp
  ParameterReader.cpp
  RandomVariable1DArray.cpp
  RandomVariable2DArray.cpp
  RandomVariable.cpp
  TableFunction.cpp
  Table.cpp
  emissionfunction.cpp
  readindata.cpp
  Poisson.cpp
  NBD.cpp
  particle_decay.cpp
  )

add_library(iSS_lib SHARED ${SOURCES})
add_executable(iSS.e ${SOURCES})
target_link_libraries (iSS_lib ${LIBS} ${ZLIB_LIBRARIES})
target_link_libraries (iSS.e ${LIBS} ${ZLIB_LIBRARIES})

install(TARGETS iSS.e DESTINATION ${CMAKE_HOME_DIRECTORY})
