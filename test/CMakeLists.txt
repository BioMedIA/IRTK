PROJECT(IRTK_TEST)

INCLUDE(CTest)

LINK_LIBRARIES(rview++ registration++ transformation++ image++ geometry++ common++ recipes)

#add google test
SET (GTEST_SOURCEDIR ../gtest)
ADD_SUBDIRECTORY(${GTEST_SOURCEDIR} ${GTEST_SOURCEDIR}/build)
INCLUDE_DIRECTORIES(${GTEST_SOURCEDIR}/include ${source_dir})

ENABLE_TESTING()

# files containing the actual tests
SET(test_sources
    applications/combineLabels_test.cc
    geometry++/irtkPointSet_test.cc
    geometry++/irtkAffineTransform_test.cc
    applications/makevolume_test.cc
)

# other source files dependencies
SET ( source_dependencies
    ${IRTK_SOURCE_DIR}/applications/src/combineLabels_core.cc
    ${IRTK_SOURCE_DIR}/applications/src/makevolume_core.cc
)

INCLUDE_DIRECTORIES ( ${IRTK_SOURCE_DIR}
                      ${IRTK_INCLUDE_DIRS}
                      ${IRTK_SOURCE_DIR}/applications/include
                    )
ADD_EXECUTABLE(irtk_tests ${test_sources} ${source_dependencies})
TARGET_LINK_LIBRARIES(irtk_tests gtest_main)
ADD_TEST(irtk_tests irtk_tests)
