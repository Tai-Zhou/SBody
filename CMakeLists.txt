cmake_minimum_required(VERSION 3.12)
project(SBody VERSION 0.0.6.220525)
add_compile_definitions(VERSION="${PROJECT_VERSION}" HAVE_INLINE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -Wshadow -Wuninitialized")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O1 -g -fsanitize=address -fno-optimize-sibling-calls")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(CMAKE_HOST_SYSTEM_NAME STREQUAL Linux)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fPIC")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
    add_compile_definitions(SBODY_RELEASE)
else()
    add_compile_definitions(BOOST_UBLAS_CHECK_RANGE)
endif()

###############
# Options
###############

option(WITH_OPENMP "Enable OpenMP" ON)
option(WITH_PYTHON "Enable Python binding" ON)
option(BUILD_DOC_ONLY "Build documentation only" OFF)

###############
# FetchContent dependencies
###############

include(FetchContent)

# CLI11
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.6.2
)

# fmt
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 11.2.0
)

# indicators
FetchContent_Declare(
    indicators
    GIT_REPOSITORY https://github.com/p-ranav/indicators.git
    GIT_TAG v2.3
)

FetchContent_MakeAvailable(CLI11 fmt indicators)

###############
# Boost
###############

find_package(Boost 1.71 REQUIRED COMPONENTS unit_test_framework)

if(NOT BOOST_FOUND)
    message(FATAL_ERROR "Please install Boost")
endif()

###############
# OpenMP
###############


if(WITH_OPENMP)
    find_package(OpenMP)

    if(OpenMP_FOUND)
        add_compile_definitions(WITH_OPENMP)
        list(APPEND OPTIONAL_LIBS OpenMP::OpenMP_CXX)
    else()
        set(WITH_OPENMP OFF)
        message(NOTICE "OpenMP not found")
    endif()
endif()

###############
# Doxygen
###############

find_package(Doxygen)

if(BUILD_DOC_ONLY AND NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen not found, please install it")
endif()

if(DOXYGEN_FOUND)
    set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/doxygen")
    set(DOXYGEN_GENERATE_HTML NO)
    set(DOXYGEN_GENERATE_LATEX NO)
    set(DOXYGEN_GENERATE_XML YES)
    set(DOXYGEN_XML_OUTPUT .)
    doxygen_add_docs(docs "${PROJECT_SOURCE_DIR}/include")
    if (BUILD_DOC_ONLY)
        return()
    endif()
endif()

###############
# Core library
###############

add_library(sbody_static STATIC
    src/file/file.cpp
    src/file/writer.cpp
    src/file/csv.cpp
    src/file/numpy.cpp
    src/infra/progress.cpp
    src/infra/status.cpp
)

target_include_directories(
    sbody_static
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(
    sbody_static
    PUBLIC
        fmt::fmt
        indicators
)

###############
# Main executable
###############

add_executable(
    sbody
    src/main.cpp
)

target_link_libraries(
    sbody
    PRIVATE
        sbody_static
        CLI11::CLI11
        Boost::boost
        ${OPTIONAL_LIBS}
)

###############
# Python binding
###############

if(WITH_PYTHON)
    add_subdirectory(python)
endif()

###############
# Tests
###############

enable_testing()

add_subdirectory(test)

add_test(
    NAME help-message
    COMMAND sbody -h
)
