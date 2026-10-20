set(TEST_SOURCES
    math/geometry_test.cpp
    math/special_functions_test.cpp
    math/solver/polynomial_test.cpp
    math/solver/system/nonlinear_system_test.cpp
    metric/inertial_frame_test.cpp
    metric/kerr_eom_test.cpp
    metric/kerr_metric_test.cpp
    metric/schwarzschild_eom_test.cpp
)

foreach(TEST_FILE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_include_directories(${TEST_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(${TEST_NAME} PUBLIC sbody_static Boost::unit_test_framework fmt::fmt)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
