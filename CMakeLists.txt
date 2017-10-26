CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/python.cmake)

SET(PROJECT_NAME parametric-curves)
SET(PROJECT_DESCRIPTION "Template based classes for creating and manipulating parametric curves. Comes with extra options specific to end-effector trajectories in robotics." )

SET(PROJECT_URL "")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/build/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib/")
# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
SET(CMAKE_VERBOSE_MAKEFILE True)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
SET(BOOST_COMPONENTS unit_test_framework serialization)
SEARCH_FOR_BOOST()

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" OFF)
IF(BUILD_PYTHON_INTERFACE)
# search for python
	FINDPYTHON(2.7 REQUIRED)
	find_package( PythonLibs 2.7 REQUIRED )
	include_directories( ${PYTHON_INCLUDE_DIRS} )
        SET(BOOST_COMPONENTS ${BOOST_COMPONENTS} python)
	#include_directories( ${Boost_INCLUDE_DIR} )

	add_subdirectory (python)

ENDIF(BUILD_PYTHON_INTERFACE)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
add_subdirectory (src/tests/spline_test)

install(FILES
#	${CMAKE_SOURCE_DIR}/include/parametriccurves/bernstein.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/curve_abc.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/splines/exact_cubic.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/MathDefs.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/splines/polynom.h
#	${CMAKE_SOURCE_DIR}/include/parametriccurves/spline_deriv_constraint.h
#	${CMAKE_SOURCE_DIR}/include/parametriccurves/bezier_curve.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/splines/cubic_spline.h
#	${CMAKE_SOURCE_DIR}/include/parametriccurves/curve_constraint.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/splines/quintic_spline.h
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametriccurves
	)
    
#install(FILES
#	${CMAKE_SOURCE_DIR}/include/parametriccurves/helpers/effector_spline.h
#	${CMAKE_SOURCE_DIR}/include/parametriccurves/helpers/effector_spline_rotation.h
#	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametriccurves/helpers
#	)

SETUP_PROJECT()
SETUP_PROJECT_FINALIZE()
