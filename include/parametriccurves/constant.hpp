/**
* \file sinusoid.hpp
* \brief Generates InfiniteSinusoidal trajectory
* \author Rohan Budhiraja
* \date 2017
*
*/

#ifndef _parameteric_curves_constant_hpp
#define _parameteric_curves_constant_hpp

#include <parametriccurves/abstract-curve.hpp>

namespace parametriccurves
{

/// \class InfiniteSinusoid
/// \brief Creates InfiniteSinusoid curve
/// The sinusoid is actually a cosine so that it starts with zero velocity.
/// Returns x = x_init + A*cos(2*pi*f*t) where f is give by 1/(2*traj_time)
template<typename Numeric=double, std::size_t Dim=1,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>>
struct Constant :
    public AbstractCurve<Numeric, Point>
{
  typedef Point          point_t;
  typedef Numeric 	 time_t;
  typedef Numeric	 num_t;
  typedef AbstractCurve<Numeric, Point> curve_abc_t;

public:
  ///\brief Constructor

  Constant()
    : curve_abc_t(-1,-1),
      x_init(point_t::Zero()) {}

  Constant(const point_t& x_init_= Eigen::Matrix<Numeric, Dim, 1>::Zero())
    : curve_abc_t(-1,-1),
      x_init(x_init_) {}

  ///\brief Destructor
  ~Constant(){}
  
public:

  virtual const point_t& operator()(const time_t& t) const
  {
    return x_init;
  }

  virtual const point_t& derivate(const time_t& t, const std::size_t& order) const
  {
    return point_t::Zero()
  }

  virtual bool setInitialPoint(const point_t& x_init_)
  {
    if(x_init.size()!=x_init_.size())
      return false;
    x_init = x_init_;
  }

  virtual bool setInitialPoint(const double& x_init_)
  {
    if(Dim!=1)
      return false;
    x_init[0] = x_init_;
    return true;
  }

protected:
  /*Attributes*/
  point_t x_init;
};
}
#endif //_CLASS_EXACTCUBIC
