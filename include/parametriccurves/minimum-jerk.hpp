/**
* \file minimum-jerk.hpp
* \brief Generates MinimumJerk trajectory
* \author Rohan Budhiraja
* \date 2017
*
*/

#ifndef _parameteric_curves_minimum_jerk_hpp
#define _parameteric_curves_minimum_jerk_hpp

#include <parametriccurves/abstract-curve.hpp>

namespace parametriccurves
{

      /** \class MinimumJerk
       *  \brief Creates MinimumJerk curve
       */
template<typename Numeric=double, std::size_t Dim=1,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>>
struct MinimumJerk :
    public AbstractCurve<Numeric, Point>
{
  typedef Point         point_t;
  typedef Numeric 	time_t;
  typedef Numeric	num_t;
  typedef AbstractCurve<Numeric, Point> curve_abc_t;

public:
  ///\brief Constructor

  MinimumJerk()
    : curve_abc_t(-1,-1),
      x_init(point_t::Zero()),
      x_final(point_t::Zero())
  {
  }

  MinimumJerk(const time_t& traj_time_,
              const point_t& x_init_= Eigen::Matrix<Numeric, Dim, 1>::Zero(),
              const point_t& x_final_= Eigen::Matrix<Numeric, Dim, 1>::Zero())
    : curve_abc_t(0,traj_time_),
      x_init(x_init_),
      x_final(x_final_)
  {
  }

  MinimumJerk(const time_t& traj_time_,
              const point_t& x_init_= Eigen::Matrix<Numeric, Dim, 1>::Zero(),
              const point_t& x_final_= Eigen::Matrix<Numeric, Dim, 1>::Zero())
    : curve_abc_t(0,traj_time_),
      x_init(x_init_),
      x_final(x_final_)  {}

  ///\brief Destructor
  ~MinimumJerk(){}
  
public:

  virtual point_t operator()(const time_t& t) const
  {
    time_t td  = t/this->t_max;
    time_t td2 = td*td;    time_t td3 = td2*td;
    time_t td4 = td3*td;   time_t td5 = td4*td;
    time_t p   = 10*td3 - 15*td4 + 6*td5;
    return x_init + (x_final-x_init)*p;
  }

  virtual point_t derivate(const time_t& t, const std::size_t& order) const
  {
    time_t td  = t/this->t_max;
    time_t td2 = td*td;    time_t td3 = td2*td;
    time_t td4 = td3*td;   time_t td5 = td4*td;
    if(order==1) {
      time_t dp  = (30*td2 - 60*td3 + 30*td4)/this->t_max;
      return (x_final-x_init)*dp;
    }
    else if(order==2) {
      time_t ddp = (60*td - 180*td2 + 120*td3)/(this->t_max*this->t_max);
      return (x_final-x_init)*ddp;
    }
    else {
      std::cerr<<"Higher order derivatives not supported"<<std::endl;
      return point_t::Zero(Dim);
    }
  }

public:
  /*Setters*/
  virtual bool setInitialPoint(const point_t& x_init_)
  {
    if(x_init_.size()!=x_init.size())
      return false;
    x_init = x_init_;
  }

  virtual bool setInitialPoint(const num_t& x_init_)
  {
    if(Dim!=1)
      return false;
    x_init[0] = x_init_;
    return true;
  }

  virtual bool setFinalPoint(const point_t& x_final_)
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

protected:
  /*Attributes*/
  point_t x_init;
  point_t x_final;

private:

};
}
#endif //_CLASS_EXACTCUBIC
