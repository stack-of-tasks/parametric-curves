cmake_minimum_required(VERSION 3.14)

# Project properties
set(PROJECT_ORG stack-of-tasks)
set(PROJECT_NAME parametric-curves)
set(PROJECT_DESCRIPTION
    "Template based classes for creating and manipulating parametric curves. Comes with extra options specific to end-effector trajectories in robotics."
)
set(PROJECT_URL "http://github.com/${PROJECT_ORG}/${PROJECT_NAME}")

# Project options
option(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)
option(INSTALL_PYTHON_INTERFACE_ONLY "Install *ONLY* the python bindings" OFF)

# Project configuration
if(NOT INSTALL_PYTHON_INTERFACE_ONLY)
  set(PROJECT_USE_CMAKE_EXPORT TRUE)
endif(NOT INSTALL_PYTHON_INTERFACE_ONLY)
set(CUSTOM_HEADER_DIR ${PROJECT_NAME})
set(CXX_DISABLE_WERROR TRUE)
set(DOXYGEN_USE_MATHJAX YES)

# Check if the submodule cmake have been initialized
set(JRL_CMAKE_MODULES "${CMAKE_CURRENT_LIST_DIR}/cmake")
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/cmake/base.cmake")
  message(STATUS "JRL cmakemodules not found. Let's fetch it.")
  include(FetchContent)
  FetchContent_Declare(
    "jrl-cmakemodules"
    GIT_REPOSITORY "https://github.com/jrl-umi3218/jrl-cmakemodules.git")
  FetchContent_MakeAvailable("jrl-cmakemodules")
  FetchContent_GetProperties("jrl-cmakemodules" SOURCE_DIR JRL_CMAKE_MODULES)
endif()

include("${JRL_CMAKE_MODULES}/base.cmake")
include("${JRL_CMAKE_MODULES}/boost.cmake")

set_default_cmake_build_type("RelWithDebInfo")

# Project definition
compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})

# Project dependencies
if(BUILD_PYTHON_INTERFACE)
  add_project_dependency(eigenpy 2.7.12 REQUIRED)
endif(BUILD_PYTHON_INTERFACE)

add_project_dependency(Eigen3 REQUIRED)
add_project_dependency(Boost REQUIRED COMPONENTS serialization)

if(BUILD_PYTHON_INTERFACE)
  string(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
  set(PYWRAP ${PY_NAME}_pywrap)
  add_subdirectory(python)
endif(BUILD_PYTHON_INTERFACE)

# Main Library
set(${PROJECT_NAME}_HEADERS
    # include/${PROJECT_NAME}/bernstein.h
    include/${PROJECT_NAME}/abstract-curve.hpp
    include/${PROJECT_NAME}/MathDefs.h
    include/${PROJECT_NAME}/spline.hpp
    include/${PROJECT_NAME}/polynomial.hpp
    include/${PROJECT_NAME}/infinite-sinusoid.hpp
    include/${PROJECT_NAME}/infinite-const-acc.hpp
    include/${PROJECT_NAME}/constant.hpp
    include/${PROJECT_NAME}/linear-chirp.hpp
    include/${PROJECT_NAME}/minimum-jerk.hpp
    include/${PROJECT_NAME}/text-file.hpp
    # include/${PROJECT_NAME}/bezier_curve.h
    include/${PROJECT_NAME}/curve-constraint.hpp
    include/${PROJECT_NAME}/serialization/eigen-matrix.hpp
    # include/${PROJECT_NAME}/serialization/archive.hpp
    include/${PROJECT_NAME}/spatial/force-curve.hpp
    include/${PROJECT_NAME}/utils/file-io.hpp)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} SYSTEM
                           INTERFACE ${EIGEN3_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME}
                           INTERFACE $<INSTALL_INTERFACE:include>)

if(NOT INSTALL_PYTHON_BINDINGS_ONLY)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${TARGETS_EXPORT_NAME}
    DESTINATION lib)
endif(NOT INSTALL_PYTHON_BINDINGS_ONLY)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif(BUILD_TESTING)
