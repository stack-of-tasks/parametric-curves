/**
* \file linear-chirp.hpp
* \brief Generates LinearChirp trajectory
* \author Rohan Budhiraja
* \date 2017
*
*/

#ifndef _parameteric_curves_linear_chirp_hpp
#define _parameteric_curves_linear_chirp_hpp

#include <parametriccurves/abstract-curve.hpp>

namespace parametriccurves
{

      /** \class LinearChirp
       *  \brief Creates LinearChirp curve
       *  Linear chirp trajectory generator.
       *  A linear chirp is a sinusoid whose frequency is a linear function of time.
       *  In particular the frequency starts from a value f0 and it increases linearly
       *  up to a value f1. Then it goes back to f0 and the trajectory is ended.
       */
template<typename Numeric=double, std::size_t Dim=1,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>>
struct LinearChirp :
    public AbstractCurve<Numeric, Point>
{
  typedef Point         point_t;
  typedef Point         freq_t;
  typedef Point         phase_t;
  typedef Numeric 	time_t;
  typedef Numeric	num_t;
  typedef AbstractCurve<Numeric, Point> curve_abc_t;

public:
  ///\brief Constructor

  LinearChirp()
    : curve_abc_t(-1,-1),
      x_init(point_t::Zero()),
      x_final(point_t::Zero())
  {
    f0.setZero(Dim);
    f1.setZero(Dim);
  }

  LinearChirp(const time_t& traj_time_,
              const point_t& x_init_= Eigen::Matrix<Numeric, Dim, 1>::Zero(),
              const point_t& x_final_= Eigen::Matrix<Numeric, Dim, 1>::Zero())
    : curve_abc_t(0,traj_time_),
      x_init(x_init_),
      x_final(x_final_)
  {
    f0.setZero(Dim);
    f1.setZero(Dim);
  }

  LinearChirp(const time_t& traj_time_,
              const freq_t& f0_,
              const freq_t& f1_,
              const point_t& x_init_= Eigen::Matrix<Numeric, Dim, 1>::Zero(),
              const point_t& x_final_= Eigen::Matrix<Numeric, Dim, 1>::Zero())
    : curve_abc_t(0,traj_time_),
      f0(f0_),
      f1(f1_),
      x_init(x_init_),
      x_final(x_final_)  {}

  ///\brief Destructor
  ~LinearChirp(){}
  
public:

  virtual point_t operator()(const time_t& t) const
  {
    const point_t& m_p   = 0.5*(1.0-phase(t).array().cos());
    return x_init.array() + (x_final.array()-x_init.array())*m_p.array();
  }

  virtual point_t derivate(const time_t& t, const std::size_t& order) const
  {
    if(order==1) {
      const point& m_dp  = M_PI*freq(t).array() * phase(t).array().sin();
      return (x_final-x_init).array()*m_dp.array();
    }
    else if(order==2) {
      const point& m_ddp = 2.0*M_PI*M_PI* freq(t).array()* freq(t).array()* phase(t).array().cos();
      return (x_final-x_init).array()*m_ddp.array();
    }
    else {
      std::cerr<<"Higher order derivatives not supported"<<std::endl;
      return point_t::Zero(Dim);
    }
  }

public:
  /*Setters*/
  virtual bool setInitialFrequency(const Eigen::VectorXd& f0)
  {
    if(f0.size()!=f0_.size())
      return false;
    f0 = f0_;
    return true;
  }

  virtual bool setInitialFrequency(const double& f0_)
  {
    if(Dim!=1)
      return false;
    f0[0] = f0_;
    return true;
  }

  virtual bool setFinalFrequency(const Eigen::VectorXd& f1_)
  {
    if(f1.size()!=f1_.size())
      return false;
    f1 = f1_;
    return true;
  }


  virtual bool setFinalFrequency(const double& f1_)
  {
    if(Dim!=1)
      return false;
    f1[0] = f1_;
    return true;
  }

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

protected:
  /*Attributes*/
  freq_t f0;          /// initial frequency
  freq_t f1;          /// final frequency
  point_t x_init;
  point_t x_final;

private:
  inline static const num_t& freq(const time_t& t) const 
  {
    const freq_t& m_k = 2.0*(f1-f0)/this->t_max;
    if(t<0.5*this->t_max)
      return f0 + m_k*t;
    else
      return f1 + m_k*(0.5*this->t_max - t);
  }

  inline static const num_t& phase(const time_t& t) const 
  {
    const freq_t& m_k = 2.0*(f1-f0)/this->t_max;
    const phase_t& m_phi_0 = M_PI*this->t_max*(f0-f1);
    if(t<0.5*this->t_max)
      return 2*M_PI*t*(f0 + 0.5*m_k*t);
    else
      return m_phi_0 + 2*M_PI*t*(f1 + 0.5*m_k*(this->t_max - t));
  }

};
}
#endif //_CLASS_EXACTCUBIC
