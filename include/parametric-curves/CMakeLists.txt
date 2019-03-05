install(FILES
  #	bernstein.h
  abstract-curve.hpp
  MathDefs.h
  spline.hpp
  polynomial.hpp
  infinite-sinusoid.hpp
  infinite-const-acc.hpp
  constant.hpp
  linear-chirp.hpp
  minimum-jerk.hpp
  text-file.hpp
  #	bezier_curve.h
  curve-constraint.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves
  )

install(FILES
  serialization/eigen-matrix.hpp
  #${CMAKE_SOURCE_DIR}/include/parametric-curves/serialization/archive.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves/serialization
  )
install(FILES
  spatial/force-curve.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves/spatial
  )
install(FILES
  utils/file-io.hpp
  #${CMAKE_SOURCE_DIR}/include/parametric-curves/serialization/archive.hpp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/parametric-curves/utils
  )
