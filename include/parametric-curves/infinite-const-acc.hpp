/**
* \file constant-acceleration.hpp
* \brief Generates InfiniteConstAcc trajectory
* \author Rohan Budhiraja
* \date 2017
*
*/

#ifndef _parameteric_curves_infinite_constant_acceleration_hpp
#define _parameteric_curves_infinite_constant_acceleration_hpp

#include <parametric-curves/abstract-curve.hpp>
#include <cmath>

namespace parametriccurves
{

/// \class InfiniteConstAcc
/// \brief Creates InfiniteConstAcc curve
/// s = s_0 + u_0*t+0.5*a_0*t^2

template<typename Numeric=double, Eigen::Index Dim=1,
         typename Point= Eigen::Matrix<Numeric, Dim, 1> >
struct InfiniteConstAcc :
    public AbstractCurve<Numeric, Point>
{
  typedef Point          point_t;
  typedef Numeric 	time_t;
  typedef Numeric	num_t;

  typedef AbstractCurve<Numeric, Point> curve_abc_t;

public:
  ///\brief Constructor

  InfiniteConstAcc(const time_t& traj_time_,
                       const point_t& x_init_= Eigen::Matrix<Numeric, Dim, 1>::Zero(),
                       const point_t& x_final_= Eigen::Matrix<Numeric, Dim, 1>::Zero())
    : curve_abc_t(-1,-1),
      traj_time(traj_time_),
      x_init(x_init_),
      x_final(x_final_) {}

  ///\brief Destructor
  ~InfiniteConstAcc(){}
  
public:

  virtual const point_t operator()(const time_t& t_) const
  {
    const time_t t = std::fmod(t_, 2*traj_time);
    const point_t& m_ddx0 = 4.0*(x_final-x_init)/(traj_time*traj_time);
    if(t<=0.5*traj_time)
      return x_init + 0.5*m_ddx0*t*t;
    else if(t>0.5*traj_time && t<=1.5*traj_time)
      return (x_init+x_final)/2
        +0.5*m_ddx0*traj_time*(t-0.5*traj_time)
        -0.5*m_ddx0*(t-0.5*traj_time)*(t-0.5*traj_time);
    else //(t>1.5*traj_time && t<=2*traj_time)
      return (x_init+x_final)/2
        - 0.5*m_ddx0*traj_time*(t-1.5*traj_time)
        + 0.5*m_ddx0*(t-1.5*traj_time)*(t-1.5*traj_time);
  }
  
  virtual const point_t derivate(const time_t& t_, const std::size_t& order) const
  {
    const time_t t = std::fmod(t_, 2*traj_time);
    const point_t& m_ddx0 = 4.0*(x_final-x_init)/(traj_time*traj_time);
    if(order==1) {
      if(t<=0.5*traj_time)
        return m_ddx0*t;
      else if(t>0.5*traj_time && t<=1.5*traj_time)
        return 0.5*m_ddx0*traj_time - m_ddx0*(t-0.5*traj_time);
      else //(t>1.5*traj_time && t<=2*traj_time)
        return -0.5*m_ddx0*traj_time + m_ddx0*(t-1.5*traj_time);
    }
    else if(order==2) {
      if(t<=0.5*traj_time)
        return m_ddx0;
      else if(t>0.5*traj_time && t<=1.5*traj_time)
        return -m_ddx0;
      else //(t>1.5*traj_time && t<=2*traj_time)
        return m_ddx0;
    }
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
};
}
#endif //_CLASS_EXACTCUBIC
