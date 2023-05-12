add_unit_test(spline_tests Main.cpp)
target_link_libraries(spline_tests ${PROJECT_NAME})

if(BUILD_PYTHON_INTERFACE)
  # add_subdirectory(python)
endif(BUILD_PYTHON_INTERFACE)
