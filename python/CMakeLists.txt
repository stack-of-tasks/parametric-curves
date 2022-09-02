add_library(${PYWRAP} SHARED spline_python.cpp)
target_link_libraries(${PYWRAP} eigenpy::eigenpy Boost::serialization)

set_target_properties(${PYWRAP} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PY_NAME})

if(APPLE)
  # We need to change the extension for python bindings
  set_target_properties(${PYWRAP} PROPERTIES SUFFIX ".so")
endif(APPLE)

install(TARGETS ${PYWRAP} DESTINATION "${PYTHON_SITELIB}/${PY_NAME}")

python_install_on_site(${PY_NAME} "__init__.py")
