CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/python.cmake)

SET(PROJECT_NAME parametriccurves)
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
  FINDPYTHON(2.7 REQUIRED)
  INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
ENDIF(BUILD_PYTHON_INTERFACE)

SET(BOOST_COMPONENTS ${BOOST_REQUIERED_COMPONENTS} ${BOOST_OPTIONAL_COMPONENTS} ${BOOST_BUILD_COMPONENTS})
SEARCH_FOR_BOOST()
SEARCH_FOR_EIGEN()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

install(FILES
#	${CMAKE_SOURCE_DIR}/include/parametriccurves/bernstein.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/abstract-curve.hpp
	${CMAKE_SOURCE_DIR}/include/parametriccurves/MathDefs.h
	${CMAKE_SOURCE_DIR}/include/parametriccurves/spline.hpp
	${CMAKE_SOURCE_DIR}/include/parametriccurves/polynomial.hpp
	${CMAKE_SOURCE_DIR}/include/parametriccurves/infinite-sinusoid.hpp
        #	${CMAKE_SOURCE_DIR}/include/parametriccurves/spline_deriv_constraint.h
        #	${CMAKE_SOURCE_DIR}/include/parametriccurves/bezier_curve.h
        ${CMAKE_SOURCE_DIR}/include/parametriccurves/curve-constraint.hpp
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametriccurves
	)
    
install(FILES
	${CMAKE_SOURCE_DIR}/include/parametriccurves/serialization/eigen-matrix.hpp
	#${CMAKE_SOURCE_DIR}/include/parametriccurves/serialization/archive.hpp
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametriccurves/serialization
        )

IF(BUILD_PYTHON_INTERFACE)
  SET(PYWRAP ${PROJECT_NAME}_pywrap)
  ADD_SUBDIRECTORY(python)
ENDIF(BUILD_PYTHON_INTERFACE)

ADD_SUBDIRECTORY (src/tests/spline_test)


PKG_CONFIG_APPEND_BOOST_LIBS(${BOOST_REQUIERED_COMPONENTS})
SETUP_PROJECT_FINALIZE()
