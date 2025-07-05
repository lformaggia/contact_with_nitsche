# clean_all.cmake

# Helper macro to empty a directory (but keep the directory itself)
macro(empty_directory dir_path)
    if(EXISTS "${dir_path}")
        file(GLOB CONTENTS "${dir_path}/*")
        if(CONTENTS)
            set(FILES_TO_REMOVE "")
            foreach(item IN LISTS CONTENTS)
                get_filename_component(filename "${item}" NAME)
                if(NOT filename STREQUAL "report.pdf")
                    list(APPEND FILES_TO_REMOVE "${item}")
                endif()
            endforeach()
            if(FILES_TO_REMOVE)
                file(REMOVE_RECURSE ${FILES_TO_REMOVE})
                message(STATUS "Emptied directory (except report.pdf): ${dir_path}")
            else()
                message(STATUS "Directory only contains report.pdf, nothing removed: ${dir_path}")
            endif()
        else()
            message(STATUS "Directory already empty: ${dir_path}")
        endif()
    else()
        message(STATUS "Directory does not exist: ${dir_path}")
    endif()
endmacro()


# Clean build lib and doc directories
empty_directory("${PROJECT_ROOT}/build")
empty_directory("${PROJECT_ROOT}/lib")
empty_directory("${PROJECT_ROOT}/doc")


# Clean each example directory
file(GLOB EXAMPLES RELATIVE "${PROJECT_ROOT}/examples" "${PROJECT_ROOT}/examples/*")

foreach(EXAMPLE ${EXAMPLES})
    set(EXAMPLE_DIR "${PROJECT_ROOT}/examples/${EXAMPLE}")

    # Skip non-directories
    if(NOT IS_DIRECTORY "${EXAMPLE_DIR}")
        continue()
    endif()

    file(GLOB EX_FILES "${EXAMPLE_DIR}/*")

    foreach(FILE ${EX_FILES})
        get_filename_component(FILE_NAME "${FILE}" NAME)

        if(NOT FILE_NAME STREQUAL "data.pot" AND NOT FILE_NAME STREQUAL "output")
            file(REMOVE_RECURSE "${FILE}")
        endif()
    endforeach()

    # Clear output folder if it exists
    if(EXISTS "${EXAMPLE_DIR}/output")
        file(GLOB OUTPUT_FILES "${EXAMPLE_DIR}/output/*")
        foreach(OUT_FILE ${OUTPUT_FILES})
            file(REMOVE_RECURSE "${OUT_FILE}")
        endforeach()
    endif()
endforeach()


message(STATUS "Running clean_all.cmake")
