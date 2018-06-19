CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/python.cmake)

SET(PROJECT_NAME parametric-curves)
SET(PROJECT_DESCRIPTION "Template based classes for creating and manipulating parametric curves. Comes with extra options specific to end-effector trajectories in robotics." )
SET(PROJECT_URL "http://github.com/stack-of-tasks/parameteric-curves")


#set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/build/")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")
#set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib/")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
SET(CMAKE_VERBOSE_MAKEFILE True)

SETUP_PROJECT()

IF(WIN32)
  SET(LINK copy_if_different)
ELSE(WIN32)
  SET(LINK create_symlink)
ENDIF(WIN32)


ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.0.5")
SET(BOOST_REQUIERED_COMPONENTS serialization)
SET(BOOST_OPTIONAL_COMPONENTS "")

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
IF(BUILD_PYTHON_INTERFACE)
  SET(BOOST_OPTIONAL_COMPONENTS ${BOOST_OPTIONAL_COMPONENTS} python)
  FINDPYTHON()
  INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
ENDIF(BUILD_PYTHON_INTERFACE)

SET(BOOST_COMPONENTS ${BOOST_REQUIERED_COMPONENTS} ${BOOST_OPTIONAL_COMPONENTS} ${BOOST_BUILD_COMPONENTS})
SEARCH_FOR_BOOST()
SEARCH_FOR_EIGEN()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

install(FILES
  #	${CMAKE_SOURCE_DIR}/include/parametric-curves/bernstein.h
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/abstract-curve.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/MathDefs.h
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/spline.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/polynomial.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/infinite-sinusoid.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/infinite-const-acc.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/constant.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/linear-chirp.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/minimum-jerk.hpp
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/text-file.hpp
  #	${CMAKE_SOURCE_DIR}/include/parametric-curves/bezier_curve.h
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/curve-constraint.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves
  )

install(FILES
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/serialization/eigen-matrix.hpp
  #${CMAKE_SOURCE_DIR}/include/parametric-curves/serialization/archive.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves/serialization
  )
install(FILES
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/spatial/force-curve.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves/spatial
  )
install(FILES
  ${CMAKE_SOURCE_DIR}/include/parametric-curves/utils/file-io.hpp
  #${CMAKE_SOURCE_DIR}/include/parametric-curves/serialization/archive.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves/utils
  )

IF(BUILD_PYTHON_INTERFACE)
  STRING(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
  SET(PYWRAP ${PY_NAME}_pywrap)
  ADD_SUBDIRECTORY(python)
ENDIF(BUILD_PYTHON_INTERFACE)

ADD_SUBDIRECTORY (src/tests/spline_test)


PKG_CONFIG_APPEND_BOOST_LIBS(${BOOST_REQUIERED_COMPONENTS})
SETUP_PROJECT_FINALIZE()
