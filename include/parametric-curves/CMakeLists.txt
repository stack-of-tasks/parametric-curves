install(FILES
  #bernstein.h
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
  #bezier_curve.h
  curve-constraint.hpp
  DESTINATION include/${CUSTOM_HEADER_DIR}
  )

install(FILES
  serialization/eigen-matrix.hpp
  #serialization/archive.hpp
  #serialization/archive.hpp
  DESTINATION include/${CUSTOM_HEADER_DIR}/serialization
  )
install(FILES
  spatial/force-curve.hpp
  DESTINATION include/${CUSTOM_HEADER_DIR}/spatial
  )
install(FILES
  utils/file-io.hpp
  DESTINATION include/${CUSTOM_HEADER_DIR}/utils
  )
