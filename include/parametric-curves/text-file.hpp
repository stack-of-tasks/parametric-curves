/**
 * \file text-file.hpp
 * \brief Reads from text file
 *
 */

#ifndef _parameteric_curves_text_file_hpp
#define _parameteric_curves_text_file_hpp

#include <parametric-curves/abstract-curve.hpp>
#include <parametric-curves/utils/file-io.hpp>

namespace parametriccurves {

/// \class TextFile.
/// \brief Loads curve from file
template <typename Numeric = double, Eigen::Index Dim = Eigen::Dynamic,
          typename Point = Eigen::Matrix<Numeric, Dim, 1> >
struct TextFile : public AbstractCurve<Numeric, Point> {
  typedef Point point_t;
  typedef Numeric time_t;
  typedef Numeric num_t;
  typedef Eigen::Matrix<Numeric, Dim, Eigen::Dynamic> pos_t;
  typedef Eigen::Matrix<Numeric, Dim, Eigen::Dynamic> vel_t;
  typedef Eigen::Matrix<Numeric, Dim, Eigen::Dynamic> acc_t;

  typedef AbstractCurve<Numeric, Point> curve_abc_t;

 public:
  ///\brief Constructor

  TextFile(const time_t& dt_, const std::size_t& size_)
      : curve_abc_t(-1, -1), timeStep(dt_), size(size_) {}

  ///\brief Destructor
  ~TextFile() {}

 public:
  virtual const point_t operator()(const time_t& t) const {
    Eigen::VectorXd::Index i = (Eigen::VectorXd::Index)std::floor(t / timeStep);
    return posValues.row(i);
  }

  virtual const point_t derivate(const time_t& t,
                                 const std::size_t& order) const {
    Eigen::VectorXd::Index i = (Eigen::VectorXd::Index)std::floor(t / timeStep);
    if (order == 1)
      return velValues.row(i);
    else if (order == 2)
      return accValues.row(i);
    else {
      std::cerr << "Higher order derivatives not supported" << std::endl;
      return point_t::Zero(size);
    }
  }

 public:
  virtual bool loadTextFile(const std::string& fileName) {
    Eigen::MatrixXd data =
        parametriccurves::utils::readMatrixFromFile(fileName);
    if (data.cols() == size) {
      std::cout << fileName << ": setting derivatives to zero" << std::endl;
      posValues = data;
      velValues.setZero(data.rows(), size);
      accValues.setZero(data.rows(), size);
    } else if (data.cols() == 2 * size) {
      std::cout << fileName << ": setting second derivative to zero"
                << std::endl;
      posValues = data.leftCols(size);
      velValues = data.rightCols(size);
      accValues = accValues.setZero(data.rows(), size);
    } else if (data.cols() == 3 * size) {
      posValues = data.leftCols(size);
      velValues = data.middleCols(size, size);
      accValues = data.rightCols(size);
    } else {
      std::cout << "Unexpected number of columns (expected " << size << " or "
                << 2 * size << " or " << 3 * size << ", found " << data.cols()
                << ")\n";
      return false;
    }
    this->t_max = timeStep * (double)data.rows();
    this->t_min = 0.0;
    x_init = posValues.row(0);
    return true;
  }

  virtual bool setInitialPoint(const point_t& /*x_init*/) { return false; }
  virtual bool setInitialPoint(const num_t& /*x_init*/) { return false; }

  const point_t& getInitialPoint(void) const { return x_init; }

 protected:
  /*Attributes*/
  point_t x_init;
  pos_t posValues;
  vel_t velValues;
  acc_t accValues;
  time_t timeStep;
  std::size_t size;
};
}  // namespace parametriccurves
#endif  //_CLASS_EXACTCUBIC
