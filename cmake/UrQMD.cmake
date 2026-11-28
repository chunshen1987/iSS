# Build UrQMD as a library and keep its source checkout untouched.
set(ISS_URQMD_SOURCE_DIR "${PROJECT_SOURCE_DIR}/urqmd_code/urqmd/src"
    CACHE PATH "Directory containing the UrQMD 3.4 Fortran sources")
enable_language(Fortran)
find_package(Threads REQUIRED)
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "The bundled UrQMD 3.4 build requires GNU Fortran")
endif()
find_program(ISS_PATCH_EXECUTABLE patch)
if(NOT ISS_PATCH_EXECUTABLE)
    message(FATAL_ERROR "The optional UrQMD build requires the patch command")
endif()

if(NOT EXISTS "${ISS_URQMD_SOURCE_DIR}/urqmd.f")
    set(urqmd_checkout "${CMAKE_CURRENT_BINARY_DIR}/_deps/urqmd_afterburner")
    set(urqmd_fetched_source_dir "${urqmd_checkout}/urqmd/src")
    if(NOT EXISTS "${urqmd_fetched_source_dir}/urqmd.f"
       AND NOT EXISTS "${urqmd_checkout}/.git")
        find_package(Git REQUIRED)
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_deps")
        message(STATUS "UrQMD sources not found in ${ISS_URQMD_SOURCE_DIR}. Cloning from Bitbucket.")
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" clone --depth 1
                    "https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git"
                    "${urqmd_checkout}"
            RESULT_VARIABLE urqmd_clone_result
            ERROR_VARIABLE urqmd_clone_error)
        if(NOT urqmd_clone_result STREQUAL "0")
            message(FATAL_ERROR
                "Could not clone UrQMD into ${urqmd_checkout}\n${urqmd_clone_error}\n"
                "Set ISS_URQMD_SOURCE_DIR to an existing UrQMD 3.4 source directory to build offline.")
        endif()
    endif()
    if(NOT EXISTS "${urqmd_fetched_source_dir}/urqmd.f")
        message(FATAL_ERROR
            "Expected UrQMD sources at ${urqmd_fetched_source_dir}. "
            "Use a fresh build directory or set ISS_URQMD_SOURCE_DIR to a valid local checkout.")
    endif()
    set(ISS_URQMD_SOURCE_DIR "${urqmd_fetched_source_dir}")
endif()

set(urqmd_generated_dir "${CMAKE_CURRENT_BINARY_DIR}/urqmd")
file(MAKE_DIRECTORY "${urqmd_generated_dir}")
foreach(unit urqmd init)
    configure_file("${ISS_URQMD_SOURCE_DIR}/${unit}.f"
                   "${urqmd_generated_dir}/${unit}.f" COPYONLY)
endforeach()
set(urqmd_patch "${CMAKE_CURRENT_LIST_DIR}/urqmd-library.patch")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${urqmd_patch}")
execute_process(COMMAND "${ISS_PATCH_EXECUTABLE}" -f -N -l -p1 -i "${urqmd_patch}"
                WORKING_DIRECTORY "${urqmd_generated_dir}"
                RESULT_VARIABLE urqmd_patch_result
                OUTPUT_VARIABLE urqmd_patch_output ERROR_VARIABLE urqmd_patch_error)
if(NOT urqmd_patch_result EQUAL 0)
    message(FATAL_ERROR "Could not adapt the UrQMD 3.4 event loop:\n${urqmd_patch_output}${urqmd_patch_error}")
endif()

# Source list from UrQMD's GNUmakefile. Files included by these units are
# resolved from the original source directory.
set(urqmd_units
    addpart.f angdis.f anndec.f bessel.f blockres.f boxprg.f cascinit.f
    coload.f dectim.f delpart.f detbal.f dwidth.f error.f getmass.f getspin.f
    hepnam.f hepchg.f hepcmp.f input.f iso.f ityp2pdg.f jdecay2.f make22.f
    numrec.f output.f paulibl.f proppot.f pythia6409.f quadri.f90 saveinfo.f
    scatter.f siglookup.f string.f tabinit.f upmerge.f whichres.f gnuranf.f)
set(urqmd_sources)
foreach(unit IN LISTS urqmd_units)
    list(APPEND urqmd_sources "${ISS_URQMD_SOURCE_DIR}/${unit}")
endforeach()
add_library(iss_urqmd_fortran STATIC ${urqmd_sources}
            "${urqmd_generated_dir}/urqmd.f" "${urqmd_generated_dir}/init.f"
            "${PROJECT_SOURCE_DIR}/src/urqmd_bridge.f")
set_target_properties(iss_urqmd_fortran PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(iss_urqmd_fortran PRIVATE "${ISS_URQMD_SOURCE_DIR}")
target_compile_options(iss_urqmd_fortran PRIVATE -cpp -std=legacy)
