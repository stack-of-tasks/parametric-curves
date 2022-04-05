/**
 * \file polynomial.hpp
 * \brief Definition of a cubic spline.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 *
 * This file contains definitions for the Polynomial struct.
 * It allows the creation and evaluation of natural
 * smooth splines of arbitrary dimension and order
 */

#ifndef _parameteric_curves_polynomial_hpp
#define _parameteric_curves_polynomial_hpp

#include <Eigen/Dense>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <parametric-curves/abstract-curve.hpp>
#include <parametric-curves/serialization/eigen-matrix.hpp>
#include <stdexcept>
#include <vector>

namespace parametriccurves {
/// \class Polynomial
/// \brief Represents a Polynomialf arbitrary order defined on the interval
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + ... + d(t - t_min_)^N, where N is the order
///
template <typename Numeric = double, Eigen::Index Dim = 3,
          typename Point = Eigen::Matrix<Numeric, Dim, 1> >
struct Polynomial : public parametriccurves::AbstractCurve<Numeric, Point> {
  typedef Point point_t;
  typedef Numeric time_t;
  typedef Numeric num_t;
  typedef std::vector<Point, Eigen::aligned_allocator<Point> > t_point_t;
  typedef AbstractCurve<Numeric, Point> curve_abc_t;
  typedef Eigen::Matrix<double, Dim, Eigen::Dynamic> coeff_t;
  typedef Eigen::Ref<coeff_t> coeff_t_ref;

 public:
  ///\brief Constructor
  ///\param coefficients : a reference to an Eigen matrix where each column is a
  /// coefficient,
  /// from the zero order coefficient, up to the highest order. Spline order is
  /// given by the number of the columns -1.
  ///\param min: LOWER bound on interval definition of the spline
  ///\param max: UPPER bound on interval definition of the spline
  Polynomial(const coeff_t& coefficients, const time_t tmin, const time_t tmax)
      : curve_abc_t(tmin, tmax),
        coefficients_(coefficients),
        dim_(Dim),
        order_(coefficients_.cols() - 1) {
    safe_check();
  }

  ///\brief Constructor
  ///\param coefficients : a container containing all coefficients of the
  /// spline, starting
  /// with the zero order coefficient, up to the highest order. Spline order is
  /// given by the size of the coefficients
  ///\param min: LOWER bound on interval definition of the spline
  ///\param max: UPPER bound on interval definition of the spline
  Polynomial(const t_point_t& coefficients, const time_t tmin,
             const time_t tmax)
      : curve_abc_t(tmin, tmax),
        coefficients_(init_coeffs(coefficients.begin(), coefficients.end())),
        dim_(Dim),
        order_(coefficients_.cols() - 1) {
    safe_check();
  }

  Polynomial() {}

  ///\brief Constructor
  ///\param zeroOrderCoefficient : an iterator pointing to the first element of
  /// a structure containing the coefficients
  /// it corresponds to the zero degree coefficient
  ///\param out   : an iterator pointing to the last element of a structure
  /// ofcoefficients \param min: LOWER bound on interval definition of the
  /// spline \param max: UPPER bound on interval definition of the spline
  template <typename In>
  Polynomial(In zeroOrderCoefficient, In out, const time_t tmin,
             const time_t tmax)
      : curve_abc_t(tmin, tmax),
        coefficients_(init_coeffs(zeroOrderCoefficient, out)),
        dim_(Dim),
        order_(coefficients_.cols() - 1) {
    safe_check();
  }

  ///\brief Destructor
  ~Polynomial() {
    // NOTHING
  }

  Polynomial(const Polynomial& other)
      : curve_abc_t(other.t_min, other.t_max),
        coefficients_(other.coefficients_),
        dim_(other.dim_),
        order_(other.order_) {}

  // Polynomial& operator=(const Polynomial& other);

 private:
  void safe_check() {
    if (this->t_min > this->t_max) std::out_of_range("TODO");
    if (static_cast<std::size_t>(coefficients_.size()) != order_ + 1)
      std::runtime_error("Spline order and coefficients do not match");
  }

  /* Constructors - destructors */

  /*Operations*/
 public:
  /*///  \brief Evaluation of the cubic spline at time t.
    ///  \param t : the time when to evaluate the spine
    ///  \param return : the value x(t)
    virtual point_t operator()(const time_t t) const
    {
        if((t < t_min_ || t > t_max_) && Safe){ throw
    std::out_of_range("TODO");} time_t const dt (t-t_min_); time_t cdt(1);
        point_t currentPoint_ = point_t::Zero();
        for(int i = 0; i < order_+1; ++i, cdt*=dt)
            currentPoint_ += cdt *coefficients_.col(i);
        return currentPoint_;
    }*/

  ///  \brief Evaluation of the cubic spline at time t using horner's scheme.
  ///  \param t : the time when to evaluate the spine
  ///  \param return : the value x(t)
  virtual const point_t operator()(const time_t& t) const {
    if ((t < this->t_min || t > this->t_max)) {
      throw std::out_of_range("TODO");
    }
    const time_t& dt(t);
    point_t h = coefficients_.col(order_);
    std::size_t i = order_ - 1;
    bool ok = true;
    if (order_ != 0) {
      while (ok && order_ != 0) {
        h = dt * h + coefficients_.col(i);
        if (i == 0)
          ok = false;
        else
          i--;
      }
      return h;
    } else
      return h;
  }

  ///  \brief Evaluation of the derivative spline at time t.
  ///  \param t : the time when to evaluate the spline
  ///  \param order : order of the derivative
  ///  \param return : the value x(t)
  virtual const point_t derivate(const time_t& t,
                                 const std::size_t& order) const {
    if ((t < this->t_min || t > this->t_max)) {
      throw std::out_of_range("TODO");
    }
    const time_t& dt(t);
    time_t cdt(1);
    point_t currentPoint_ = point_t::Zero(dim_);
    for (std::size_t i = order; i < order_ + 1; ++i, cdt *= dt)
      currentPoint_ += cdt * coefficients_.col(i) * fact(i, order);
    return currentPoint_;
  }

  virtual const std::size_t& size() const { return dim_; }

  virtual bool setInitialPoint(const point_t& /*x_init*/) { return false; }
  virtual bool setInitialPoint(const num_t& /*x_init*/) { return false; }

 protected:
  coeff_t coefficients_;  // const
  std::size_t dim_;       // const
  std::size_t order_;     // const

 private:
  // Serialization of the class
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive& ar, const unsigned int /*version*/) const {
    ar& boost::serialization::make_nvp("dim", dim_);
    ar& boost::serialization::make_nvp("order", order_);
    ar& boost::serialization::make_nvp("coefficients", coefficients_);
    ar& boost::serialization::make_nvp("t_min", this->t_min);
    ar& boost::serialization::make_nvp("t_max", this->t_max);
  }

  template <class Archive>
  void load(Archive& ar, const unsigned int /*version*/) {
    ar& boost::serialization::make_nvp("dim", dim_);
    ar& boost::serialization::make_nvp("order", order_);
    ar& boost::serialization::make_nvp("coefficients", coefficients_);
    ar& boost::serialization::make_nvp("t_min", this->t_min);
    ar& boost::serialization::make_nvp("t_max", this->t_max);
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  num_t fact(const std::size_t n, const std::size_t order) const {
    std::size_t res(1);
    for (std::size_t i = 0; i < order; ++i) res *= n - i;
    return static_cast<num_t>(res);
  }

  /*Attributes*/
  template <typename In>
  coeff_t init_coeffs(In zeroOrderCoefficient, In highestOrderCoefficient) {
    std::size_t size =
        std::distance(zeroOrderCoefficient, highestOrderCoefficient);
    coeff_t res = coeff_t(Dim, size);
    int i = 0;
    for (In cit = zeroOrderCoefficient; cit != highestOrderCoefficient;
         ++cit, ++i)
      res.col(i) = *cit;
    return res;
  }
};  // class Polynomial
}  // namespace parametriccurves
#endif  //_STRUCT_POLYNOMIAL
