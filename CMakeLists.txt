cmake_minimum_required(VERSION 3.12)
project(SBody VERSION 0.0.6.220525)
add_compile_definitions(VERSION="${PROJECT_VERSION}" HAVE_INLINE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(CMAKE_HOST_SYSTEM_NAME STREQUAL Linux)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_definitions(GSL_RANGE_CHECK_OFF)
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(SBODY_LIB_SHARED OFF)
set(WITH_CFITSIO OFF)
set(WITH_OPENMP ON)
set(WITH_PYTHON OFF)
set(RECORD_TRACE OFF)
set(VIEW_TAU OFF)
set(VIEW_HAMILTONIAN ON)

find_package(GSL REQUIRED)

if(NOT GSL_FOUND)
    message(FATAL_ERROR "Please install GSL")
endif()

find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/include fmt/include indicators/include ${Python3_INCLUDE_DIRS})
add_subdirectory(${PROJECT_SOURCE_DIR}/fmt)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_executable(sbody src/SBody.cpp)
target_link_libraries(sbody SBodyLibStatic)

find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/doxygen")
    set(DOXYGEN_GENERATE_HTML NO)
    set(DOXYGEN_GENERATE_LATEX NO)
    set(DOXYGEN_GENERATE_XML YES)
    set(DOXYGEN_XML_OUTPUT .)
    doxygen_add_docs(docs "${PROJECT_SOURCE_DIR}/include")
endif()
