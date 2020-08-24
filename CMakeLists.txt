CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# Project properties
SET(PROJECT_ORG stack-of-tasks)
SET(PROJECT_NAME parametric-curves)
SET(PROJECT_DESCRIPTION "Template based classes for creating and manipulating parametric curves. Comes with extra options specific to end-effector trajectories in robotics." )
SET(PROJECT_URL "http://github.com/${PROJECT_ORG}/${PROJECT_NAME}")

# Project options
OPTION(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)
OPTION(INSTALL_PYTHON_INTERFACE_ONLY "Install *ONLY* the python bindings" OFF)

# Project configuration
IF(NOT INSTALL_PYTHON_INTERFACE_ONLY)
  SET(PROJECT_USE_CMAKE_EXPORT TRUE)
ENDIF(NOT INSTALL_PYTHON_INTERFACE_ONLY)
SET(CUSTOM_HEADER_DIR ${PROJECT_NAME})
SET(CXX_DISABLE_WERROR TRUE)
SET(DOXYGEN_USE_MATHJAX YES)

# JRL-cmakemodule setup
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/python.cmake)

# Project definition
COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})

# Project dependencies
ADD_PROJECT_DEPENDENCY(Eigen3 REQUIRED)
ADD_PROJECT_DEPENDENCY(Boost REQUIRED COMPONENTS serialization)

IF(BUILD_PYTHON_INTERFACE)
  FINDPYTHON()
  SEARCH_FOR_BOOST_PYTHON()
  ADD_PROJECT_DEPENDENCY(eigenpy REQUIRED)
ENDIF(BUILD_PYTHON_INTERFACE)

IF(BUILD_PYTHON_INTERFACE)
  STRING(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
  SET(PYWRAP ${PY_NAME}_pywrap)
  ADD_SUBDIRECTORY(python)
ENDIF(BUILD_PYTHON_INTERFACE)

# Main Library
SET(${PROJECT_NAME}_HEADERS
  #include/${PROJECT_NAME}/bernstein.h
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
  #include/${PROJECT_NAME}/bezier_curve.h
  include/${PROJECT_NAME}/curve-constraint.hpp
  include/${PROJECT_NAME}/serialization/eigen-matrix.hpp
  #include/${PROJECT_NAME}/serialization/archive.hpp
  #include/${PROJECT_NAME}/serialization/archive.hpp
  include/${PROJECT_NAME}/spatial/force-curve.hpp
  include/${PROJECT_NAME}/utils/file-io.hpp
  )

ADD_LIBRARY(${PROJECT_NAME} INTERFACE)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} INTERFACE $<INSTALL_INTERFACE:include>)

IF(NOT INSTALL_PYTHON_BINDINGS_ONLY)
  INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)
ENDIF(NOT INSTALL_PYTHON_BINDINGS_ONLY)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING)
