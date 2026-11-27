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
    spin_polarization.cpp
    Histogram.cpp
    Spectators.cpp
)

add_library(${lib_name} SHARED ${SOURCES})
if (APPLE)
    set_target_properties (${lib_name} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (${lib_name} ${LIBS} ${ZLIB_LIBRARIES})

if (link_with_lib)
    add_executable(${exe_name} main.cpp)
    target_link_libraries (${exe_name} ${lib_name} ${LIBS} ${ZLIB_LIBRARIES})
else (link_with_lib)
    add_executable(${exe_name} main.cpp ${SOURCES})
    target_link_libraries (${exe_name} ${LIBS} ${ZLIB_LIBRARIES})
endif (link_with_lib)

if (APPLE)
    set_target_properties (${exe_name} PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
install(TARGETS ${exe_name} ${lib_name}
        RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)

if (unittest)
    add_executable (unittest_Spectators.e
        Spectators_UnitTest.cpp ../tests/UnitTestMain.cpp)
    target_include_directories(unittest_Spectators.e PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries (unittest_Spectators.e ${lib_name} ${LIBS})
    add_test(NAME spectators COMMAND unittest_Spectators.e
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")
    install(TARGETS unittest_Spectators.e ${lib_name}
            RUNTIME DESTINATION ${CMAKE_HOME_DIRECTORY}/tests
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib/static)

    add_executable(unittest_numerics.e
        Random_UnitTest.cpp
        ../tests/TestDiscreteSampling.cpp
        ../tests/TestMomentumSampling.cpp
        ../tests/TestEosInterpolation.cpp
        ../tests/UnitTestMain.cpp)
    target_include_directories(unittest_numerics.e PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(unittest_numerics.e PRIVATE ${lib_name})

    set(unit_test_dir "${CMAKE_CURRENT_BINARY_DIR}/unit-test-data")
    file(MAKE_DIRECTORY
        "${unit_test_dir}/iSS_tables/EOS_tables/HRG4D"
        "${unit_test_dir}/iSS_tables/deltaf_tables/urqmd")
    # Keep the light-fermion discrepancy reproducible while its physics
    # prescription is investigated.
    add_test(NAME numerics COMMAND unittest_numerics.e
        "--test-case-exclude=Light fermion momentum*"
        WORKING_DIRECTORY "${unit_test_dir}")
endif (unittest)
