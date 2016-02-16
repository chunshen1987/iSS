add_executable (iSS.e
  main.cpp 
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
  )
target_link_libraries (iSS.e ${LIBS})

install(TARGETS iSS.e DESTINATION ${CMAKE_HOME_DIRECTORY})
