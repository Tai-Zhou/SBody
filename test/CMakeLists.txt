set(TEST_SOURCES
    TestUtility.cpp
)

foreach(TEST_FILE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_NAME} SBodyLibStatic Boost::unit_test_framework)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
