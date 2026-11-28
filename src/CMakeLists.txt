set (SOURCES
    iSS.cpp
    arsenal.cpp
    ParameterReader.cpp
    RandomVariable1DArray.cpp
    RandomVariable2DArray.cpp
    Random.cpp
    TableFunction.cpp
    Table.cpp
    eos_4D.cpp
    emissionfunction.cpp
    readindata.cpp
    particle_decay.cpp
    pretty_ostream.cpp
    MomentumSamplerShell.cpp
    BosonMomentumSampler.cpp
    FermionMomentumSampler.cpp
    MomentumSamplerBase.cpp
    FSSW.cpp
    urqmd_output.cpp
    urqmd_transport.cpp
    spin_polarization.cpp
    Histogram.cpp
    Spectators.cpp
)

add_library(${lib_name} SHARED ${SOURCES})
if (APPLE)
    set_target_properties (${lib_name} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (${lib_name} ${LIBS} ${ZLIB_LIBRARIES})
if (ISS_ENABLE_URQMD)
    target_compile_definitions(${lib_name} PRIVATE ISS_WITH_URQMD)
    target_link_libraries(${lib_name} iss_urqmd_fortran Threads::Threads)
endif ()

if (link_with_lib)
    add_executable(${exe_name} main.cpp)
    target_link_libraries (${exe_name} ${lib_name} ${LIBS} ${ZLIB_LIBRARIES})
else (link_with_lib)
    add_executable(${exe_name} main.cpp ${SOURCES})
    target_link_libraries (${exe_name} ${LIBS} ${ZLIB_LIBRARIES})
    if (ISS_ENABLE_URQMD)
        target_compile_definitions(${exe_name} PRIVATE ISS_WITH_URQMD)
        target_link_libraries(${exe_name} iss_urqmd_fortran Threads::Threads)
    endif ()
endif (link_with_lib)

if (APPLE)
    set_target_properties (${exe_name} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
install(TARGETS ${exe_name} ${lib_name}
        RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)

if (unittest)
    add_executable (unittest_Spectators.e Spectators_UnitTest.cpp)
    target_link_libraries (unittest_Spectators.e ${lib_name} ${LIBS})
    install(TARGETS unittest_Spectators.e ${lib_name}
            RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}/tests
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib/static)
endif (unittest)

if (ISS_ENABLE_URQMD)
    add_executable(test_urqmd_transport ../tests/TestUrQMDTransport.cpp)
    target_include_directories(test_urqmd_transport PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(test_urqmd_transport ${lib_name})
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/urqmd-test")
    add_test(NAME urqmd_transport
             COMMAND test_urqmd_transport
             WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/urqmd-test")
    set_tests_properties(urqmd_transport PROPERTIES ENVIRONMENT "URQMD_TAB=tables.dat")

    find_program(ISS_PYTHON_EXECUTABLE python3)
    if (ISS_PYTHON_EXECUTABLE)
        add_test(NAME urqmd_in_memory_cli
                 COMMAND "${ISS_PYTHON_EXECUTABLE}"
                         "${PROJECT_SOURCE_DIR}/tests/TestUrQMDTransport.py"
                         $<TARGET_FILE:${exe_name}>)
    endif ()
endif ()
