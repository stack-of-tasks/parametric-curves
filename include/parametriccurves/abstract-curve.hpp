/**
* \file AbstractCurve.hpp
* \brief interface for a Curve of arbitrary dimension.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* Interface for a curve
*/
#include <cstddef>

#ifndef _parameteric_curves_abstract_curve_hpp
#define _parameteric_curves_abstract_curve_hpp

namespace parametriccurves
{
/// \struct AbstractCurve
/// \brief Represents a curve of dimension Dim
/// is Safe is false, no verification is made on the evaluation of the curve.
template<typename Numeric, typename Point >
struct  AbstractCurve
{
  typedef Point   point_t;
  typedef Numeric time_t;
  typedef Numeric num_t;
public:
  /* Constructors - destructors */
  AbstractCurve(time_t t_min_, time_t t_max_):  t_min(t_min_), t_max(t_max_) { }
  AbstractCurve(){ }
  virtual ~AbstractCurve(){}  
public:
  
  ///  \brief Evaluation of the cubic spline at time t.
  ///  \param t : the time when to evaluate the spine
  ///  \param return : the value x(t)
  virtual const point_t operator()(const time_t& t) const = 0;

  ///  \brief Evaluation of the derivative spline at time t.
  ///  \param t : the time when to evaluate the spline
  ///  \param order : order of the derivative
  ///  \param return : the value x(t)
  virtual const point_t derivate(const time_t& t, const std::size_t& order) const = 0;

public:
  /*Getters*/
  virtual const time_t tmin() const { return t_min; }
  virtual const time_t tmax() const { return t_max; }
  virtual bool checkRange(const time_t t) const { return (t>=t_min)&&(t<=t_max); }

  /* Setters */
  virtual bool setInitialPoint(const point_t& /*x_init*/) = 0;
  virtual bool setInitialPoint(const num_t& /*x_init*/) = 0;

  virtual bool setTimePeriod(const time_t& traj_time_)
  {
    t_min = 0.0;
    t_max = traj_time_;  
    return true;
  }

protected:
  time_t t_min;
  time_t t_max;

};
}
#endif //_STRUCT_CURVE_ABC
