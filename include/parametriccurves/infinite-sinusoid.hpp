/**
* \file sinusoid.hpp
* \brief Generates InfiniteSinusoidal trajectory
* \author Rohan Budhiraja
* \date 2017
*
*/

#ifndef _parameteric_curves_sinusoid_hpp
#define _parameteric_curves_sinusoid_hpp

#include <parametriccurves/abstract-curve.hpp>

namespace parametriccurves
{

/// \class InfiniteSinusoid
/// \brief Creates InfiniteSinusoid curve
/// The sinusoid is actually a cosine so that it starts with zero velocity.
/// Returns x = x_init + A*cos(2*pi*f*t) where f is give by 1/(2*traj_time)
template<typename Numeric=double, std::size_t Dim=1,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>>
struct InfiniteSinusoid :
    public AbstractCurve<Numeric, Point>
{
  typedef Point          point_t;
  typedef Numeric 	time_t;
  typedef Numeric	num_t;

  typedef AbstractCurve<Numeric, Point> curve_abc_t;

public:
  ///\brief Constructor

  InfiniteSinusoid()
    : curve_abc_t(-1,-1),
      traj_time(NAN),
      x_init(point_t::Zero()),
      x_final(point_t::Zero()) {}

  InfiniteSinusoid(const time_t& traj_time_,
                   const point_t& x_init_= Eigen::Matrix<Numeric, Dim, 1>::Zero(),
                   const point_t& x_final_= Eigen::Matrix<Numeric, Dim, 1>::Zero())
    : curve_abc_t(-1,-1),
      traj_time(traj_time_),
      x_init(x_init_),
      x_final(x_final_) {}

  ///\brief Destructor
  ~InfiniteSinusoid(){}
  
public:

  virtual point_t operator()(const time_t& t) const
  {
    return x_init + 0.5*(x_final-x_init)*(1.0-cos(two_pi_f(t)));
  }

  virtual point_t derivate(const time_t& t, const std::size_t& order) const
  {
    if(order==1)
      return (x_final-x_init)*0.5*two_pi_f(1)*sin(two_pi_f(t));
    else if(order==2)
      return (x_final-x_init)*0.5*two_pi_f(1)*two_pi_f(1)*cos(two_pi_f(t));
    else {
      std::cerr<<"Higher order derivatives not supported"<<std::endl;
      return point_t::Zero(Dim);
    }
  }

public:
  virtual bool setInitialPoint(const point_t& x_init_)
  {
    if(x_init_.size()!=x_init.size())
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

  virtual bool setFinalPoint(const Eigen::VectorXd& x_final_)
  {
    if(x_final.size()!=x_final_.size())
      return false;
    x_final = x_final_;
    return true;
  }

  virtual bool setFinalPoint(const double& x_final_)
  {
    if(Dim!=1)
      return false;
    x_final[0] = x_final_;
    return true;
  }
  virtual bool setTrajectoryTime(double traj_time_)
  {
    if(traj_time_<=0.0)
      return false;
    traj_time = traj_time_;
    return true;
  }

protected:
  /*Attributes*/
  point_t x_init;
  point_t x_final;
  time_t traj_time;

private:
  inline static const num_t& two_pi_f(const time_t& t) const 
  {
    return (M_PI/m_traj_time)*t;
  }

};
}
#endif //_CLASS_EXACTCUBIC
