cmake_minimum_required(VERSION 3.16)
project(PACSProject CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(PACS_PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Output folders
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# -------------------------------
# Options
# -------------------------------
option(ENABLE_EXAMPLES "Enable running example tests via CTest" OFF)
option(EXAMPLES_VERBOSE "Pass -v to examples" OFF)
option(EXAMPLES_USE_GMSH "Pass -m to examples" OFF)

# -------------------------------
# Step 1: Build and install muparserx
# -------------------------------
include(ExternalProject)

set(MUPARSERX_INSTALL_DIR ${CMAKE_SOURCE_DIR}/lib/muparserx)
set(MUPARSERX_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/muparserx)
set(MUPARSERX_BINARY_DIR ${CMAKE_BINARY_DIR}/muparserx)

ExternalProject_Add(muparserx_ext
    SOURCE_DIR ${MUPARSERX_SOURCE_DIR}
    BINARY_DIR ${MUPARSERX_BINARY_DIR}
    INSTALL_DIR ${MUPARSERX_INSTALL_DIR}
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DBUILD_SHARED_LIBS=ON
        -DCMAKE_BUILD_TYPE=Release
    BUILD_BYPRODUCTS ${MUPARSERX_INSTALL_DIR}/lib/libmuparserx.so
)

# -------------------------------
# Step 2: Static library from src/
# -------------------------------
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${MUPARSERX_INSTALL_DIR}/include/muparserx
    /usr/local/include
    /usr/local/include/getfem
    /usr/local/include/gmm
)

file(GLOB_RECURSE CORE_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)

add_library(mycontactlib STATIC ${CORE_SOURCES})
set_target_properties(mycontactlib PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
)

add_dependencies(mycontactlib muparserx_ext)

# -------------------------------
# Step 3: Build the main executable
# -------------------------------
add_executable(main main.cpp)
target_link_libraries(main PRIVATE mycontactlib)

target_link_libraries(main PRIVATE
    ${MUPARSERX_INSTALL_DIR}/lib/libmuparserx.so
    getfem
    dmumps_seq mumps_common_seq smumps_seq cmumps_seq
    lapack blas qhull
)

target_link_options(main PRIVATE -Wl,-rpath,'$ORIGIN/../lib/muparserx/lib')

# -------------------------------
# Step 4: Copy or symlink main into each example dir
# -------------------------------
file(GLOB EXAMPLE_DIRS RELATIVE ${CMAKE_SOURCE_DIR}/examples ${CMAKE_SOURCE_DIR}/examples/*)
foreach(EXAMPLE ${EXAMPLE_DIRS})
    if(IS_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/${EXAMPLE}")
        if(UNIX)
            add_custom_command(TARGET main POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E create_symlink
                        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/main
                        ${CMAKE_SOURCE_DIR}/examples/${EXAMPLE}/main
                COMMENT "Symlinking main to ${EXAMPLE}"
            )
        else()
            add_custom_command(TARGET main POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/main
                        ${CMAKE_SOURCE_DIR}/examples/${EXAMPLE}/main
                COMMENT "Copying main to ${EXAMPLE}"
            )
        endif()
    endif()
endforeach()

# -------------------------------
# Step 5: CTest-based example execution (optional)
# -------------------------------
include(CTest)
enable_testing()

# Compose optional runtime arguments
set(RUNTIME_ARGS "")
if(EXAMPLES_VERBOSE)
    list(APPEND RUNTIME_ARGS -v)
endif()
if(EXAMPLES_USE_GMSH)
    list(APPEND RUNTIME_ARGS -m)
endif()

foreach(EXAMPLE ${EXAMPLE_DIRS})
    if(EXAMPLE STREQUAL "user")
        continue()
    endif()
    set(EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/examples/${EXAMPLE}")
    set(DATA_POT "${EXAMPLE_DIR}/data.pot")
    if(IS_DIRECTORY ${EXAMPLE_DIR} AND EXISTS ${DATA_POT})
        add_test(NAME run_${EXAMPLE}
            COMMAND ${CMAKE_COMMAND} -E env
                ${EXAMPLE_DIR}/main ${RUNTIME_ARGS}
            WORKING_DIRECTORY ${EXAMPLE_DIR}
        )

        set_tests_properties(run_${EXAMPLE} PROPERTIES
            PASS_REGULAR_EXPRESSION ".*"
            FAIL_REGULAR_EXPRESSION "Cannot open file"
        )
    endif()
endforeach()

# -------------------------------
# Step 6: Custom clean target
# -------------------------------
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_ROOT=${PACS_PROJECT_ROOT}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/clean_all.cmake
    COMMENT "Performing deep clean: build/, lib/, example outputs..."
)

# ----------------------------------------
# Step 7: Doxygen documentation generation
# ----------------------------------------
find_package(Doxygen REQUIRED)

# set(DOXYGEN_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc")
set(DOXYGEN_CONFIG_FILE "${CMAKE_BINARY_DIR}/Doxyfile")

configure_file(${CMAKE_SOURCE_DIR}/Doxyfile ${DOXYGEN_CONFIG_FILE} @ONLY)

add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG_FILE}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
)