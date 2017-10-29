/**
* \file exact_cubic.h
* \brief class allowing to create an Exact cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the CubicSpline class.
* Given a set of waypoints (x_i*) and timestep (t_i), it provides the unique set of
* cubic splines fulfulling those 4 restrictions :
* - x_i(t_i) = x_i* ; this means that the curve passes through each waypoint
* - x_i(t_i+1) = x_i+1* ;
* - its derivative is continous at t_i+1
* - its 2nd derivative is continous at t_i+1
* more details in paper "Task-Space Trajectories via Cubic Spline Optimization"
* By J. Zico Kolter and Andrew Y.ng (ICRA 2009)
*/


#ifndef _parameteric_curves_splines_cubic_splines_hpp
#define _parameteric_curves_splines_cubic_splines_hpp

#include "../abstract-curve.hpp"
//#include "quintic_spline.h"
#include "polynomial.hpp"
#include "../MathDefs.h"

#include <functional>
#include <fstream>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
namespace parametriccurves
{
namespace splines
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
Polynomial<Numeric,Dim,Point,T_Point> create_cubic(Point const& a, Point const& b,
                                                   Point const& c, Point const &d,
                                                   const Numeric min, const Numeric max)
{
    T_Point coeffs = make_cubic_vector<Point, T_Point>(a,b,c,d);
    return Polynomial<Numeric,Dim,Point,T_Point>(coeffs.begin(),coeffs.end(), min, max);
}

/// \class CubicSpline
/// \brief Represents a set of cubic splines defining a continuous function 
/// crossing each of the waypoint given in its initialization
///
template<typename Numeric=double, std::size_t Dim=3,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>,
         typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> >,
         typename SplineBase=Polynomial<Numeric, Dim, Point, T_Point> >
struct CubicSpline :
    public AbstractCurve<Numeric, Point>
{
  typedef Point          point_t;
  typedef T_Point        t_point_t;
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
  typedef Numeric 	time_t;
  typedef Numeric	num_t;
  typedef SplineBase spline_t;
  typedef typename std::vector<spline_t> t_spline_t;
  typedef typename t_spline_t::iterator it_spline_t;
  typedef typename t_spline_t::const_iterator cit_spline_t;
  typedef AbstractCurve<Numeric, Point> curve_abc_t;

public:
  ///\brief Constructor
  ///\param wayPointsBegin : an iterator pointing to the first element of a waypoint container
  ///\param wayPointsEns   : an iterator pointing to the end           of a waypoint container
  /* Constructors - destructors */
  template<typename In>
  CubicSpline(In wayPointsBegin, In wayPointsEnd)
    : curve_abc_t(),
      subSplines_(computeWayPoints<In>(wayPointsBegin, wayPointsEnd))
  {
    this->t_min = subSplines_.front().tmin();
    this->t_max = subSplines_.back().tmax();
  }
  ///\brief Constructor
  ///\param subSplines: vector of subsplines
  CubicSpline(const t_spline_t& subSplines)
    : curve_abc_t(subSplines.front().tmin(), subSplines.back().tmax()),
      subSplines_(subSplines) {}

  ///\brief Copy Constructor
  CubicSpline(const CubicSpline& other)
    : curve_abc_t(other.subSplines_.front().tmin(), other.subSplines_.front().tmax()),
      subSplines_(other.subSplines_)  {}
  
  ///\brief Destructor
  ~CubicSpline(){}
  
private:
  template<typename In>
  t_spline_t computeWayPoints(In wayPointsBegin, In wayPointsEnd) const
  {
    std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
    if(size < 1) {
      throw; // TODO
    }
    t_spline_t subSplines; subSplines.reserve(size);
    
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
      subSplines.push_back(create_cubic<Numeric,Dim,Point,T_Point>(a.row(i), b.row(i),
                                                                   c.row(i), d.row(i),
                                                                   (*it).first, (*next).first));
    }
    subSplines.push_back(create_cubic<Numeric,Dim,Point,T_Point>(a.row(size-1), b.row(size-1),
                                                                 c.row(size-1), d.row(size-1),
                                                                 (*it).first, (*it).first));
    return subSplines;
  }

public:
  virtual point_t operator()(const time_t t) const
  {
    if((t < subSplines_.front().tmin() || t > subSplines_.back().tmax())) {
      throw std::out_of_range("TODO");
    }
    for(cit_spline_t it = subSplines_.begin(); it != subSplines_.end(); ++ it) {
      if(t >= (it->tmin()) && t <= (it->tmax())) {
        return it->operator()(t);
      }
    }
  }

  virtual point_t derivate(const time_t t, const std::size_t order) const
  {
    if((t < subSplines_.front().tmin() || t > subSplines_.back().tmax())) {
      throw std::out_of_range("TODO");
    }
    for(cit_spline_t it = subSplines_.begin(); it != subSplines_.end(); ++ it) {
      if(t >= (it->tmin()) && t <= (it->tmax())) {
        return it->derivate(t, order);
      }
    }
  }

  /*Getters*/
  virtual const time_t tmin() const {return subSplines_.front().tmin();}
  virtual const time_t tmax() const {return subSplines_.back().tmax();}

protected:
  /*Attributes*/
  t_spline_t subSplines_; // const

private:

  // Serialization of the class
  friend class boost::serialization::access;  
  template<class Archive>
  void save(Archive & ar, const unsigned int /*version*/) const
  {
    ar & subSplines_;
    /*    for(typename t_spline_t::iterator it = subSplines_.begin();
        it != subSplines_.end(); ++it)
        ar & boost::serialization::make_nvp("subSplines",*it);  */
    return;
  }
  
  template<class Archive>
  void load(Archive & ar, const unsigned int /*version*/)
  {
    ar & subSplines_;

    /*    for(typename t_spline_t::iterator it = subSplines_.begin();
        it != subSplines_.end(); ++it)
        ar & boost::serialization::make_nvp("subSplines",*it); */
    this->t_min = subSplines_.front().tmin();
    this->t_max = subSplines_.back().tmax();
    return;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
  void loadSpline(const std::string & filename) throw (std::invalid_argument)
  {
    std::ifstream ifs(filename.c_str());
    if(ifs) {
      boost::archive::text_iarchive ia(ifs);
      CubicSpline cubic_spline = *static_cast<CubicSpline*>(this);
      ia >> cubic_spline;
    }
    else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }
      
  /// \brief Saved a Derived object as a text file.
  void saveSpline(const std::string & filename) const throw (std::invalid_argument)
  {
    std::ofstream ofs(filename.c_str());
    if(ofs) {
      boost::archive::text_oarchive oa(ofs);
      oa << *static_cast<const CubicSpline*>(this);
    }
    else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  //BOOST_SERIALIZATION_SPLIT_MEMBER()

};
}
}
#endif //_CLASS_EXACTCUBIC
