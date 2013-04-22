# Macro that tests and returns whether a C++ feature is present in the
# current compiler

set(CXX_CHECK_FEATURE_MODULE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(CXX_FEATURES_SUPPORTED "")
set(CXX_FEATURES_UNSUPPORTED "")

macro(CXX_PERFORM_TEST TEST_SOURCE_FILE TEST_TEST_BINARY_DIR EXPECTED_RESULT RESULT COMPILE_DEFINITIONS)

    try_run(
        RUN_RESULT_VAR COMPILE_RESULT_VAR
        "${TEST_BINARY_DIR}" "${TEST_SOURCE_FILE}"
        COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS}"
        COMPILE_OUTPUT_VARIABLE COMPILE_OUT
        RUN_OUTPUT_VARIABLE RUN_OUT
    )

    set(RESULT_VAR FALSE)

    if (COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR)
        set(RESULT_VAR TRUE)
    endif (COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR)

    if (NOT ("${RESULT_VAR}" STREQUAL "${EXPECTED_RESULT}"))
        # message ("Got ${RESULT_VAR} as a result, but ${EXPECTED_RESULT} expected")

        if (NOT ${COMPILE_RESULT_VAR})
            # message("------ compilation output ------")
            # message("${COMPILE_OUT}")
        endif (NOT ${COMPILE_RESULT_VAR})

        if (${RUN_RESULT_VAR})
            # message("---------- run output ----------")
            # message("${RUN_OUT}")
            # message("Process returned: ${RUN_RESULT_VAR}")
        endif (${RUN_RESULT_VAR})

        # message("--------------------------------")

        set (${RESULT} FALSE)

    else ()
        set (${RESULT} TRUE)

    endif ()



endmacro(CXX_PERFORM_TEST TEST_SOURCE EXPECTED_RESULT RESULT)


macro(CXX_CHECK_FEATURE CXX_VERSION FEATURE_NAME FEATURE_NUMBER RESULT_VAR COMPILE_DEFINITIONS)

    # Testing whether we have previously set the variable
    if(NOT DEFINED ${RESULT_VAR})

        set(TEST_BINARY_DIR
            "${CMAKE_CURRENT_BINARY_DIR}/cxx_${FEATURE_NUMBER}"
            )

        set(TEST_SOURCE_BASE
            "${CXX_CHECK_FEATURE_MODULE_DIR}/${CXX_VERSION}-test-${FEATURE_NAME}-${FEATURE_NUMBER}"
            )

        set(TEST_SOURCE_FILE "${TEST_SOURCE_BASE}.cpp")
        set(FAILTEST_SOURCE_FILE "${TEST_SOURCE_BASE}-fail.cpp")

        set(FEATURE_NAME
            "'${FEATURE_NAME}' (${CXX_VERSION} N${FEATURE_NUMBER})"
            )

        message(STATUS "Checking C++ support for ${FEATURE_NAME}")

        string (COMPARE EQUAL "${CMAKE_CXX_COMPILER_ID}" "Clang" CMAKE_COMPILER_IS_CLANG)
        string (COMPARE EQUAL "${CMAKE_CXX_COMPILER_ID}" "GNU" CMAKE_COMPILER_IS_GCC)

        set (ADD_COMPILE_DEFINITIONS "")

        if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG OR CMAKE_COMPILER_IS_GCC)
            set (ADD_COMPILE_DEFINITIONS "-std=c++0x")
        endif ()

        if (EXISTS ${TEST_SOURCE_FILE})
            CXX_PERFORM_TEST(${TEST_SOURCE_FILE} ${TEST_BINARY_DIR} TRUE ${RESULT_VAR} "${COMPILE_DEFINITIONS} ${ADD_COMPILE_DEFINITIONS}")
        endif ()

        if (${RESULT_VAR} AND EXISTS ${FAILTEST_SOURCE_FILE})
            CXX_PERFORM_TEST(${FAILTEST_SOURCE_FILE} ${TEST_BINARY_DIR} FALSE ${RESULT_VAR} "${COMPILE_DEFINITIONS} ${ADD_COMPILE_DEFINITIONS}")
        endif ()

        if (${RESULT_VAR})
            message(STATUS "Checking C++ support for ${FEATURE_NAME} -- works")
            set (CXX_FEATURES_SUPPORTED
                "${CXX_FEATURES_SUPPORTED} ${FEATURE_NAME} (${FEATURE_NUMBER}),"
                )

        else ()
            message(STATUS "Checking C++ support for ${FEATURE_NAME} -- not supported")
            set (CXX_FEATURES_UNSUPPORTED
                "${CXX_FEATURES_UNSUPPORTED} ${FEATURE_NAME} (${FEATURE_NUMBER}),"
                )

        endif ()

        # This would break the feature reporting on second call of cmake
        # TODO: Fix?
        # set(${RESULT_VAR} ${${RESULT_VAR}} CACHE INTERNAL "C++ support for ${FEATURE_NAME}")

    endif(NOT DEFINED ${RESULT_VAR})

endmacro(CXX_CHECK_FEATURE)

