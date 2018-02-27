/**
* \file exact_cubic.h
* \brief class allowing to create an Exact cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the Spline class.
*/


#ifndef _parameteric_curves_spline_hpp
#define _parameteric_curves_spline_hpp

#include <parametric-curves/abstract-curve.hpp>
#include <parametric-curves/curve-constraint.hpp>
#include <parametric-curves/polynomial.hpp>
#include <parametric-curves/MathDefs.h>

#include <fstream>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
namespace parametriccurves
{

/// \brief Creates coefficient vector of a cubic spline defined on the interval
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + c(t - t_min_)^2 + d(t - t_min_)^3
///
template<typename Point, typename T_Point>
T_Point make_cubic_vector(Point const& a, Point const& b, Point const& c, Point const &d)
{
    T_Point res;
    res.push_back(a);res.push_back(b);res.push_back(c);res.push_back(d);
    return res;
}

template<typename Numeric, std::size_t Dim, typename Point, typename T_Point>
Polynomial<Numeric,Dim,Point> create_cubic(Point const& a, Point const& b,
                                                   Point const& c, Point const &d,
                                                   const Numeric min, const Numeric max)
{
  T_Point coeffs = make_cubic_vector<Point, T_Point>(a,b,c,d);
    return Polynomial<Numeric,Dim,Point>(coeffs.begin(),coeffs.end(), min, max);
}
/// \brief Creates coefficient vector of a quintic spline defined on the interval
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + c(t - t_min_)^2 + d(t - t_min_)^3 + e(t - t_min_)^4  + f(t - t_min_)^5
///
template<typename Point, typename T_Point>
T_Point make_quintic_vector(Point const& a, Point const& b, Point const& c,
                   Point const &d, Point const& e, Point const& f)
{
    T_Point res;
    res.push_back(a);res.push_back(b);res.push_back(c);
    res.push_back(d);res.push_back(e);res.push_back(f);
    return res;
}

template<typename Numeric, std::size_t Dim, typename Point, typename T_Point>
Polynomial<Numeric,Dim,Point> create_quintic(Point const& a, Point const& b,
                                                     Point const& c, Point const &d,
                                                     Point const &e, Point const &f,
                                                     const Numeric min, const Numeric max)  
{
  T_Point coeffs = make_quintic_vector<Point, T_Point>(a,b,c,d,e,f);
    return Polynomial<Numeric,Dim,Point>(coeffs.begin(),coeffs.end(), min, max);
}

/// \class Spline
/// \brief Represents a set of cubic splines defining a continuous function 
/// crossing each of the waypoint given in its initialization
///
template<typename Numeric=double, std::size_t Dim=Eigen::Dynamic,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>,
         typename SplineBase=Polynomial<Numeric, Dim, Point> >
struct Spline :
    public AbstractCurve<Numeric, Point>
{
  typedef Point          point_t;
  typedef std::vector<Point,Eigen::aligned_allocator<Point> > t_point_t;
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
  typedef Numeric 	time_t;
  typedef Numeric	num_t;
  typedef SplineBase spline_t;
  typedef typename std::vector<spline_t,Eigen::aligned_allocator<spline_t> > t_spline_t;
  typedef typename t_spline_t::iterator it_spline_t;
  typedef typename t_spline_t::const_iterator cit_spline_t;
  typedef curve_constraints<point_t> spline_constraints;
  typedef AbstractCurve<Numeric, Point> curve_abc_t;

public:
  ///\brief Constructor

  Spline()
    : curve_abc_t()
  {
  }

  ///\brief Constructor
  ///\param subSplines: vector of subsplines
  Spline(const t_spline_t &subSplines)
    : curve_abc_t(subSplines.front().tmin(), subSplines.back().tmax()),
      subSplines_(subSplines) {}

  ///\brief Copy Constructor
  Spline(const Spline& other)
    : curve_abc_t(other.subSplines_.front().tmin(), other.subSplines_.front().tmax()),
      subSplines_(other.subSplines_)  {}
  
  ///\brief Destructor
  ~Spline(){}
  
public:


  /* Given a set of waypoints (x_i*) and timestep (t_i), it provides the unique set of
   * cubic splines fulfulling those 4 restrictions :
   * - x_i(t_i) = x_i* ; this means that the curve passes through each waypoint
   * - x_i(t_i+1) = x_i+1* ;
   * - its derivative is continous at t_i+1
   * - its 2nd derivative is continous at t_i+1
   * more details in paper "Task-Space Trajectories via Cubic Spline Optimization"
   * By J. Zico Kolter and Andrew Y.ng (ICRA 2009) */
  template<typename In>
  void createSplineFromWayPoints(In wayPointsBegin, In wayPointsEnd)
  {
    std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
    if(size < 1) {
      throw; // TODO
    }
    subSplines_.clear();
    subSplines_.reserve(size);
    
    // refer to the paper to understand all this.
    MatrixX h1 = MatrixX::Zero(size, size);
    MatrixX h2 = MatrixX::Zero(size, size);
    MatrixX h3 = MatrixX::Zero(size, size);
    MatrixX h4 = MatrixX::Zero(size, size);
    MatrixX h5 = MatrixX::Zero(size, size);
    MatrixX h6 = MatrixX::Zero(size, size);
    
    MatrixX a =  MatrixX::Zero(size, Dim);
    MatrixX b =  MatrixX::Zero(size, Dim);
    MatrixX c =  MatrixX::Zero(size, Dim);
    MatrixX d =  MatrixX::Zero(size, Dim);
    MatrixX x =  MatrixX::Zero(size, Dim);
    In it(wayPointsBegin), next(wayPointsBegin);
    ++next;

    for(std::size_t i(0); next != wayPointsEnd; ++next, ++it, ++i) {
      num_t const dTi((*next).first  - (*it).first);
      num_t const dTi_sqr(dTi * dTi);
      num_t const dTi_cube(dTi_sqr * dTi);
      // filling matrices values
      h3(i,i)   = -3 / dTi_sqr;
      h3(i,i+1) =  3 / dTi_sqr;
      h4(i,i)   = -2 / dTi;
      h4(i,i+1) = -1 / dTi;
      h5(i,i)   =  2 / dTi_cube;
      h5(i,i+1) = -2 / dTi_cube;
      h6(i,i)   =  1 / dTi_sqr;
      h6(i,i+1) =  1 / dTi_sqr;
      if( i+2 < size) {
        In it2(next); ++ it2;
        num_t const dTi_1((*it2).first - (*next).first);
        num_t const dTi_1sqr(dTi_1 * dTi_1);
        // this can be optimized but let's focus on clarity as long as not needed
        h1(i+1, i)   =  2 / dTi;
        h1(i+1, i+1) =  4 / dTi + 4 / dTi_1;
        h1(i+1, i+2) =  2 / dTi_1;
        h2(i+1, i)   = -6 / dTi_sqr;
        h2(i+1, i+1) = (6 / dTi_1sqr) - (6 / dTi_sqr);
        h2(i+1, i+2) =  6 / dTi_1sqr;
      }
      x.row(i)= (*it).second.transpose();
    }
    // adding last x
    x.row(size-1)= (*it).second.transpose();
    a= x;
    parametriccurves::PseudoInverse(h1);
    b = h1 * h2 * x; //h1 * b = h2 * x => b = (h1)^-1 * h2 * x
    c = h3 * x + h4 * b;
    d = h5 * x + h6 * b;
    it= wayPointsBegin, next=wayPointsBegin; ++ next;

    for(int i=0; next != wayPointsEnd; ++i, ++it, ++next) {
      Numeric min = (*it).first;
      Numeric max = (*next).first;
      Point a_ = a.row(i)-b.row(i)*min+c.row(i)*min*min-d.row(i)*min*min*min;
      Point b_ = b.row(i)-2*c.row(i)*min + 3*d.row(i)*min*min;
      Point c_ = c.row(i)-3*d.row(i)*min;
      Point d_ = d.row(i);
      subSplines_.push_back(create_cubic<Numeric,Dim,Point,t_point_t>(a_,b_,c_,d_, min, max));
    }
    Numeric min = (*it).first;
    Point a_ = a.row(size-1)-b.row(size-1)*min+c.row(size-1)*min*min-d.row(size-1)*min*min*min;
    Point b_ = b.row(size-1)-2*c.row(size-1)*min + 3*d.row(size-1)*min*min;
    Point c_ = c.row(size-1)-3*d.row(size-1)*min;
    Point d_ = d.row(size-1);
    subSplines_.push_back(create_cubic<Numeric,Dim,Point,t_point_t>(a_,b_,c_,d_,  min, min));

    this->t_min = subSplines_.front().tmin();
    this->t_max = subSplines_.back().tmax();
    return;
  }

  template<typename In>
  void createSplineFromWayPointsConstr(In wayPointsBegin, In wayPointsEnd,
                                       const spline_constraints& constraints)
  {
    std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
    if(size < 1) throw; // TODO
    subSplines_.clear();
    subSplines_.reserve(size);
    spline_constraints cons = constraints;
    In it(wayPointsBegin), next(wayPointsBegin), end(wayPointsEnd-1);
    ++next;
    for(std::size_t i(0); next != end; ++next, ++it, ++i)
      compute_one_spline<In>(it, next, cons, subSplines_);
    compute_end_spline<In>(it, next,cons, subSplines_);
    return;
  }


public:
  virtual const point_t operator()(const time_t& t) const
  {
    if((t < subSplines_.front().tmin() || t > subSplines_.back().tmax())) {
      throw std::out_of_range("t is out of range");
    }
    for(cit_spline_t it = subSplines_.begin(); it != subSplines_.end(); ++ it) {
      if(t >= (it->tmin()) && t <= (it->tmax())) {
        return it->operator()(t);
      }
    }
    const point_t dummy;
    return dummy;
  }

  virtual const point_t derivate(const time_t& t, const std::size_t& order) const
  {
    if((t < subSplines_.front().tmin() || t > subSplines_.back().tmax())) {
      throw std::out_of_range("derivative call out of range");
    }
    for(cit_spline_t it = subSplines_.begin(); it != subSplines_.end(); ++ it) {
      if(t >= (it->tmin()) && t <= (it->tmax())) {
        return it->derivate(t, order);
      }
    }

    const point_t dummy;
    return dummy;
  }

  virtual const std::size_t& size() const
  {
    return subSplines_[0].size();
  }
  const t_spline_t& getSubsplines() const
  {
    return subSplines_;
  }

  virtual bool setInitialPoint(const point_t& /*x_init*/)
  {
    return false;
  }
  virtual bool setInitialPoint(const num_t& /*x_init*/)
  {
    return false;
  }
  
protected:
  /*Attributes*/
  t_spline_t subSplines_; // const

private:

    template<typename In>
    void compute_one_spline(In wayPointsBegin, In wayPointsNext, spline_constraints& constraints, t_spline_t& subSplines) const
    {
        const point_t& a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
        const point_t& b0 = constraints.init_vel , c0 = constraints.init_acc / 2.;
        const num_t& init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
        const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt;
        const point_t d0 = (a1 - a0 - b0 * dt - c0 * dt_2) / dt_3;

        Point a_ = a0-b0*init_t+c0*init_t*init_t-d0*init_t*init_t*init_t;
        Point b_ = b0-2*c0*init_t + 3*d0*init_t*init_t;
        Point c_ = c0-3*d0*init_t;
        Point d_ = d0;
        subSplines.push_back(create_cubic<Numeric,Dim,Point,t_point_t>(a_,b_,c_,d_,
                                                                       init_t, end_t));

        constraints.init_vel = subSplines.back().derivate(end_t, 1);
        constraints.init_acc = subSplines.back().derivate(end_t, 2);
    }

    template<typename In>
    void compute_end_spline(In wayPointsBegin, In wayPointsNext, spline_constraints& constraints, t_spline_t& subSplines) const
    {
        const point_t& a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
        const point_t& b0 = constraints.init_vel, b1 = constraints.end_vel,
                       c0 = constraints.init_acc / 2., c1 = constraints.end_acc;
        const num_t& init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
        const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt, dt_4 = dt_3 * dt, dt_5 = dt_4 * dt;
        //solving a system of four linear eq with 4 unknows: d0, e0
        const point_t alpha_0 = a1 - a0 - b0 *dt -    c0 * dt_2;
        const point_t alpha_1 = b1 -      b0     - 2 *c0 * dt;
        const point_t alpha_2 = c1 -               2 *c0;
        const num_t x_d_0 = dt_3, x_d_1 = 3*dt_2, x_d_2 = 6*dt;
        const num_t x_e_0 = dt_4, x_e_1 = 4*dt_3, x_e_2 = 12*dt_2;
        const num_t x_f_0 = dt_5, x_f_1 = 5*dt_4, x_f_2 = 20*dt_3;

        point_t d, e, f;
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(3,Dim);
        rhs.row(0) = alpha_0; rhs.row(1) = alpha_1; rhs.row(2) = alpha_2;
        Eigen::Matrix3d eq  = Eigen::Matrix3d::Zero();
        eq(0,0) = x_d_0; eq(0,1) = x_e_0; eq(0,2) = x_f_0;
        eq(1,0) = x_d_1; eq(1,1) = x_e_1; eq(1,2) = x_f_1;
        eq(2,0) = x_d_2; eq(2,1) = x_e_2; eq(2,2) = x_f_2;
        rhs = eq.inverse().eval() * rhs;
        d = rhs.row(0); e = rhs.row(1); f = rhs.row(2);
        num_t min = init_t;
        Numeric min2 = min*min;
        Numeric min3 = min2*min;
        Numeric min4 = min3*min;
        Numeric min5 = min4*min;
        Point a_ = a0-b0*min+c0*min2-d*min3+e*min4-f*min5;
        Point b_ = b0-2*c0*min+3*d*min2-4*e*min3+5*f*min4;
        Point c_ = c0-3*d*min+6*e*min2-10*f*min3;
        Point d_ = d-4*e*min+10*f*min2;
        Point e_ = e-5*f*min;
        Point f_ = f;

        subSplines.push_back(create_quintic<Numeric,Dim,Point,t_point_t>
                             (a_,b_,c_,d_,e_,f_, init_t, end_t));
    }


  // Serialization of the class
  friend class boost::serialization::access;  
  template<class Archive>
  void save(Archive & ar, const unsigned int /*version*/) const
  {
    ar & subSplines_;

    return;
  }
  
  template<class Archive>
  void load(Archive & ar, const unsigned int /*version*/)
  {
    ar & subSplines_;

    this->t_min = subSplines_.front().tmin();
    this->t_max = subSplines_.back().tmax();
    return;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
  bool loadFromFile(const std::string & filename) throw (std::invalid_argument)
  {
    std::ifstream ifs(filename.c_str());
    if(ifs) {
      boost::archive::text_iarchive ia(ifs);
      Spline& cubic_spline = *static_cast<Spline*>(this);
      ia >> cubic_spline;
    }
    else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
      return false;
    }
    return true;
  }
      
  /// \brief Saved a Derived object as a text file.
  bool saveToFile(const std::string & filename) const throw (std::invalid_argument)
  {
    std::ofstream ofs(filename.c_str());
    if(ofs) {
      boost::archive::text_oarchive oa(ofs);
      oa << *static_cast<const Spline*>(this);
    }
    else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
      return false;
    }
    return true;
  }

  //BOOST_SERIALIZATION_SPLIT_MEMBER()

};
}
#endif //_CLASS_EXACTCUBIC
