/**
 * \file exact_cubic.h
 * \brief class allowing to create an Exact cubic spline.
 * \author Rohan Budhiraja
 * \version 0.1
 * \date 09/11/2017
 *
 * This file contains representation of spatial vectors as splines
 */

#ifndef _parameteric_curves_force_curve_hpp
#define _parameteric_curves_force_curve_hpp

#include <parametric-curves/abstract-curve.hpp>
#include <parametric-curves/spline.hpp>

#include <fstream>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
namespace parametriccurves {
namespace spatial {

using parametriccurves::Spline;
/// \class ForceCurve
/// \brief Representation of a spatial vector curve in the form of splines
/// Returns Plucker coordinates in the form of (Linear(3), Angular(3))
///
template <typename Numeric = double>
struct ForceCurve : public AbstractCurve<Numeric, Eigen::Matrix<Numeric, 6, 1> > {
  static const std::size_t Dim = 3;
  typedef Eigen::Matrix<Numeric, 2 * Dim, 1> force_t;
  typedef Eigen::Matrix<Numeric, 2 * Dim, 1> motion_t;
  typedef Spline<Numeric, Dim, Eigen::Matrix<Numeric, Dim, 1> > spline_lin_t;
  typedef Spline<Numeric, Dim, Eigen::Matrix<Numeric, Dim, 1> > spline_ang_t;
  typedef Numeric time_t;
  typedef Numeric num_t;
  typedef AbstractCurve<num_t, force_t> curve_abc_t;

 public:
  ///\brief Constructor

  ForceCurve() : curve_abc_t() {}

  ///\brief Constructor
  ///\param subSplines: vector of subsplines
  ForceCurve(const spline_lin_t& linPart_, const spline_ang_t& angPart_)
      : curve_abc_t(linPart_.tmin(), linPart_.tmax()),
        linPart(linPart_),
        angPart(angPart_),
        motionVector(motion_t::Zero()) {}

  ///\brief Copy Constructor
  ForceCurve(const ForceCurve& other)
      : curve_abc_t(other.linPart.tmin(), other.linPart.tmax()),
        linPart(other.linPart),
        angPart(other.angPart),
        motionVector(motion_t::Zero()) {}

  ///\brief Destructor
  ~ForceCurve() {}

 public:
  virtual const force_t operator()(const time_t& t) const {
    force_t s = force_t::Zero();
    s << linPart(t), angPart(t);
    return s;
  }

  virtual const force_t derivate(const time_t&, const std::size_t&) const {
    // TODO: Implement derivative
    return force_t::Zero(Dim);
  }

  virtual const std::size_t& size() const { return linPart.size(); }

  void setMotionVector(const motion_t& motionVector_) {
    motionVector = motionVector_;
    return;
  }

  virtual bool setInitialPoint(const force_t& /*x_init*/) { return false; }
  virtual bool setInitialPoint(const num_t& /*x_init*/) { return false; }

 protected:
  /*Attributes*/
  spline_lin_t linPart;   // const
  spline_ang_t angPart;   // const
  motion_t motionVector;  // const

 private:
  // Serialization of the class
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive& ar, const unsigned int /*version*/) const {
    ar& linPart;
    ar& angPart;

    return;
  }

  template <class Archive>
  void load(Archive& ar, const unsigned int /*version*/) {
    ar& linPart;
    ar& angPart;

    motionVector = motion_t::Zero();
    this->t_min = linPart.tmin();
    this->t_max = linPart.tmax();

    assert(this->t_min == angPart.tmin());
    assert(this->t_max == angPart.tmax());
    return;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

 public:
  bool loadFromFile(const std::string& filename) {
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      boost::archive::text_iarchive ia(ifs);
      ForceCurve& force_curve = *static_cast<ForceCurve*>(this);
      ia >> force_curve;
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
      return false;
    }
    return true;
  }

  /// \brief Saved a Derived object as a text file.
  bool saveToFile(const std::string& filename) const {
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::text_oarchive oa(ofs);
      oa << *static_cast<const ForceCurve*>(this);
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
      return false;
    }
    return true;
  }

  // BOOST_SERIALIZATION_SPLIT_MEMBER()
};
}  // namespace spatial
}  // namespace parametriccurves
#endif  //_CLASS_EXACTCUBIC
