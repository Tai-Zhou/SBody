cmake_minimum_required(VERSION 3.12)
project(SBody VERSION 0.0.6.220525)
add_compile_definitions(VERSION="${PROJECT_VERSION}" HAVE_INLINE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(SBODY_LIB_SHARED OFF)
set(WITH_CFITSIO OFF)
set(WITH_OPENMP ON)
set(WITH_PYTHON ON)
set(RECORD_TRACE OFF)
set(VIEW_TAU OFF)
set(VIEW_HAMILTONIAN ON)

find_package(GSL REQUIRED)

if(NOT GSL_FOUND)
    message(FATAL_ERROR "Please install GSL")
endif()

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/include indicators/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_executable(${PROJECT_NAME} src/SBody.cpp)
target_link_libraries(${PROJECT_NAME} SBodyLibStatic)

find_package(Doxygen REQUIRED)

set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/doxygen")
set(DOXYGEN_GENERATE_HTML NO)
set(DOXYGEN_GENERATE_LATEX NO)
set(DOXYGEN_GENERATE_XML YES)
set(DOXYGEN_XML_OUTPUT .)
doxygen_add_docs(docs "${PROJECT_SOURCE_DIR}/src")
