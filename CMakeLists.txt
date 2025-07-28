cmake_minimum_required(VERSION 3.12)
project(SBody VERSION 0.0.6.220525)
add_compile_definitions(VERSION="${PROJECT_VERSION}" HAVE_INLINE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(CMAKE_HOST_SYSTEM_NAME STREQUAL Linux)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fPIC")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_definitions(GSL_RANGE_CHECK_OFF)
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(SBODY_LIB_SHARED OFF)
set(WITH_OPENMP ON)
set(WITH_PYTHON ON)
option(BUILD_DOC_ONLY "Build documentation without compiling the program" OFF)

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

find_package(GSL REQUIRED)
find_package(Boost 1.71 REQUIRED COMPONENTS unit_test_framework)

if(NOT GSL_FOUND)
    message(FATAL_ERROR "Please install GSL")
endif()

if(NOT BOOST_FOUND)
    message(FATAL_ERROR "Please install Boost")
endif()

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

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/include fmt/include indicators/include ${Python3_INCLUDE_DIRS})
add_subdirectory(${PROJECT_SOURCE_DIR}/fmt)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/test)
add_executable(sbody src/SBody.cpp)
target_link_libraries(sbody SBodyLibStatic)

enable_testing()

add_test(help-message ${EXECUTABLE_OUTPUT_PATH}/sbody -h)
