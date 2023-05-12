// #include "parametriccurves/bezier_curve.h"
#include <parametric-curves/polynomial.hpp>
#include <parametric-curves/spatial/force-curve.hpp>
#include <parametric-curves/spline.hpp>

// #include "parametriccurves/splines/spline_deriv_constraint.h"
#include <boost/python.hpp>
#include <eigenpy/eigenpy.hpp>
#include <eigenpy/memory.hpp>
#include <parametric-curves/curve-constraint.hpp>
#include <vector>

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
namespace bp = boost::python;

typedef double real;
typedef Eigen::Vector3d point_t;
typedef Eigen::Matrix<double, 6, 1, 0, 6, 1> point6_t;
typedef Eigen::Matrix<double, 3, 1, 0, 3, 1> ret_point_t;
typedef Eigen::Matrix<double, 6, 1, 0, 6, 1> ret_point6_t;
typedef Eigen::VectorXd time_waypoints_t;
typedef Eigen::VectorXd time_vector_t;
typedef Eigen::Matrix<real, 3, Eigen::Dynamic> point_list_t;
typedef Eigen::Matrix<real, 6, Eigen::Dynamic> point_list6_t;
typedef std::vector<point_t, Eigen::aligned_allocator<point_t> > t_point_t;
typedef std::vector<point6_t, Eigen::aligned_allocator<point6_t> > t_point6_t;
typedef std::pair<real, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;
typedef std::pair<real, point6_t> Waypoint6;
typedef std::vector<Waypoint6> T_Waypoint6;
typedef std::vector<t_point_t, Eigen::aligned_allocator<t_point_t> >
    t3d_poly_coeffs_vector_t;
typedef typename t3d_poly_coeffs_vector_t::iterator it3d_poly_coeffs_vector_t;
typedef typename t3d_poly_coeffs_vector_t::const_iterator
    cit3d_poly_coeffs_vector_t;

// typedef spline::bezier_curve  <real, real, 3, true, point_t> bezier_t;
// typedef spline::bezier_curve  <real, real, 6, true, point6_t> bezier6_t;
typedef parametriccurves::Polynomial<real, 3, point_t> polynom_t;
typedef typename std::vector<polynom_t, Eigen::aligned_allocator<polynom_t> >
    t_spline_t;
typedef parametriccurves::Spline<real, 3, point_t> spline_t;
typedef parametriccurves::spatial::ForceCurve<real> force_t;
typedef polynom_t::coeff_t coeff_t;
typedef std::pair<real, point_t> waypoint_t;
typedef std::vector<waypoint_t, Eigen::aligned_allocator<point_t> >
    t_waypoint_t;

// typedef spline::spline_deriv_constraint  <real, real, 3, true, point_t,
// t_point_t> spline_deriv_constraint_t;
typedef parametriccurves::curve_constraints<point_t> curve_constraints_t;
typedef parametriccurves::curve_constraints<point6_t> curve_constraints6_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/

// EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_t)
// EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier6_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(polynom_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_t)

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curve_constraints_t)
// EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_deriv_constraint_t)

namespace parametriccurves {
using namespace boost::python;
template <typename PointList, typename T_Point>
T_Point vectorFromEigenArray(const PointList& array) {
  T_Point res;
  for (int i = 0; i < array.cols(); ++i) res.push_back(array.col(i));
  return res;
}
/*
template <typename Bezier, typename PointList, typename T_Point>
Bezier* wrapBezierConstructorTemplate(const PointList& array, const real lb =
0., const real ub =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), lb, ub);
}

template <typename Bezier, typename PointList, typename T_Point, typename
CurveConstraints> Bezier* wrapBezierConstructorConstraintsTemplate(const
PointList& array, const CurveConstraints& constraints, const real lb = 0., const
real ub =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), constraints, lb, ub);
}
*/
/*3D constructors */
/*
bezier_t* wrapBezierConstructor(const point_list_t& array)
{
    return wrapBezierConstructorTemplate<bezier_t, point_list_t,
t_point_t>(array) ;
}
bezier_t* wrapBezierConstructorBounds(const point_list_t& array, const real lb,
const real ub)
{
    return wrapBezierConstructorTemplate<bezier_t, point_list_t,
t_point_t>(array, lb, ub) ;
}
bezier_t* wrapBezierConstructorConstraints(const point_list_t& array, const
curve_constraints_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier_t, point_list_t,
t_point_t, curve_constraints_t>(array, constraints) ;
}
bezier_t* wrapBezierConstructorBoundsConstraints(const point_list_t& array,
const curve_constraints_t& constraints, const real lb, const real ub)
{
    return wrapBezierConstructorConstraintsTemplate<bezier_t, point_list_t,
t_point_t, curve_constraints_t>(array, constraints, lb, ub) ;
}
*/
/*END 3D constructors */
/*6D constructors */
/*
bezier6_t* wrapBezierConstructor6(const point_list6_t& array)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t,
t_point6_t>(array) ;
}
bezier6_t* wrapBezierConstructorBounds6(const point_list6_t& array, const real
lb, const real ub)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t,
t_point6_t>(array, lb, ub) ;
}
bezier6_t* wrapBezierConstructor6Constraints(const point_list6_t& array, const
curve_constraints6_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t,
t_point6_t, curve_constraints6_t>(array, constraints) ;
}
bezier6_t* wrapBezierConstructorBounds6Constraints(const point_list6_t& array,
const curve_constraints6_t& constraints, const real lb, const real ub)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t,
t_point6_t, curve_constraints6_t>(array, constraints, lb, ub) ;
}
*/
/*END 6D constructors */

polynom_t* wrapSplineConstructor(const coeff_t& array) {
  return new polynom_t(array, 0., 1.);
}

t_waypoint_t getWayPoints(const coeff_t& array,
                          const time_waypoints_t& time_wp) {
  t_waypoint_t res;
  for (int i = 0; i < array.cols(); ++i)
    res.push_back(std::make_pair(time_wp(i), array.col(i)));
  return res;
}

template <typename BezierType, int dim>
Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> wayPointsToList(
    const BezierType& self) {
  typedef typename BezierType::t_point_t t_point;
  typedef typename BezierType::t_point_t::const_iterator cit_point;
  const t_point& wps = self.waypoints();
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> res(dim, wps.size());
  int col = 0;
  for (cit_point cit = wps.begin(); cit != wps.end(); ++cit, ++col)
    res.block<dim, 1>(0, col) = *cit;
  return res;
}
void spline_from_waypoints(spline_t& self, const coeff_t& array,
                           const time_waypoints_t& time_wp) {
  t_waypoint_t wps = getWayPoints(array, time_wp);
  self.createSplineFromWayPoints(wps.begin(), wps.end());
  return;
}

spline_t* spline_by_concatenation_constructor(const bp::list& list_splines) {
  t_spline_t subSplines;
  subSplines.clear();
  for (int i = 0; i < len(list_splines); ++i) {
    spline_t _sp = bp::extract<spline_t>(list_splines[i]);
    const t_spline_t& _vec_subspline = _sp.getSubsplines();
    subSplines.insert(subSplines.end(), _vec_subspline.begin(),
                      _vec_subspline.end());
  }
  return new spline_t(subSplines);
}

spline_t* wrapExactCubicConstructorvoid() { return new spline_t(); }

spline_t* wrapExactCubicConstructorPolySequence(
    const bp::list& list_polynomials, const time_vector_t& time_vector) {
  typedef std::vector<t_point_t, Eigen::aligned_allocator<t_point_t> >
      t3d_poly_coeffs_vector_t;
  t3d_poly_coeffs_vector_t poly_coeffs_vector;
  t_spline_t subSplines;
  subSplines.clear();

  assert(time_vector.size() == len(list_polynomials) + 1);
  for (int i = 0; i < len(list_polynomials); ++i) {
    subSplines.push_back(polynom_t(bp::extract<coeff_t>(list_polynomials[i]),
                                   time_vector[i], time_vector[i + 1]));
    // time_vector[i], time_vector[i+1]));
  }
  return new spline_t(subSplines);
}

force_t* wrapForceCurveConstructorvoid() { return new force_t(); }

force_t* wrapForceCurveConstructorSplines(const spline_t& linear_part,
                                          const spline_t& ang_part) {
  return new force_t(linear_part, ang_part);
}

void spline_from_waypoints_constr(spline_t& self, const coeff_t& array,
                                  const time_waypoints_t& time_wp,
                                  const curve_constraints_t& constraints) {
  t_waypoint_t wps = getWayPoints(array, time_wp);
  self.createSplineFromWayPointsConstr(wps.begin(), wps.end(), constraints);
  return;
}

point_t get_init_vel(const curve_constraints_t& c) { return c.init_vel; }

point_t get_init_acc(const curve_constraints_t& c) { return c.init_acc; }

point_t get_end_vel(const curve_constraints_t& c) { return c.end_vel; }

point_t get_end_acc(const curve_constraints_t& c) { return c.end_acc; }

void set_init_vel(curve_constraints_t& c, const point_t& val) {
  c.init_vel = val;
}

void set_init_acc(curve_constraints_t& c, const point_t& val) {
  c.init_acc = val;
}

void set_end_vel(curve_constraints_t& c, const point_t& val) {
  c.end_vel = val;
}

void set_end_acc(curve_constraints_t& c, const point_t& val) {
  c.end_acc = val;
}

BOOST_PYTHON_MODULE(libparametric_curves_pywrap) {
  /** BEGIN eigenpy init**/
  bp::import("eigenpy");

  eigenpy::enableEigenPySpecific<point_t>();
  eigenpy::enableEigenPySpecific<ret_point_t>();
  eigenpy::enableEigenPySpecific<point_list_t>();
  eigenpy::enableEigenPySpecific<point6_t>();
  eigenpy::enableEigenPySpecific<ret_point6_t>();
  eigenpy::enableEigenPySpecific<point_list6_t>();
  eigenpy::enableEigenPySpecific<coeff_t>();
  /*eigenpy::exposeAngleAxis();
  eigenpy::exposeQuaternion();*/
  /** END eigenpy init**/

  /** BEGIN bezier curve 6**/
  /*
  class_<bezier6_t>
      ("bezier6", no_init)
          .def("__init__", make_constructor(&wrapBezierConstructor6))
          .def("__init__", make_constructor(&wrapBezierConstructorBounds6))
          //.def("__init__",
  make_constructor(&wrapBezierConstructor6Constraints))
          //.def("__init__",
  make_constructor(&wrapBezierConstructorBounds6Constraints)) .def("min",
  &bezier6_t::min) .def("max", &bezier6_t::max) .def("__call__",
  &bezier6_t::operator()) .def("derivate", &bezier6_t::derivate)
          .def("compute_derivate", &bezier6_t::compute_derivate)
          .def("compute_primitive", &bezier6_t::compute_primitive)
          .def("waypoints", &wayPointsToList<bezier6_t,6>)
          .def_readonly("degree", &bezier6_t::degree_)
          .def_readonly("nbWaypoints", &bezier6_t::size_)
      ;
  */
  /** END bezier curve**/

  /** BEGIN bezier curve**/
  /*
  class_<bezier_t>
      ("bezier", no_init)
          .def("__init__", make_constructor(&wrapBezierConstructor))
          .def("__init__", make_constructor(&wrapBezierConstructorBounds))
          .def("__init__", make_constructor(&wrapBezierConstructorConstraints))
          .def("__init__",
  make_constructor(&wrapBezierConstructorBoundsConstraints)) .def("min",
  &bezier_t::min) .def("max", &bezier_t::max) .def("__call__",
  &bezier_t::operator()) .def("derivate", &bezier_t::derivate)
          .def("compute_derivate", &bezier_t::compute_derivate)
          .def("compute_primitive", &bezier_t::compute_primitive)
          .def("waypoints", &wayPointsToList<bezier_t,3>)
          .def_readonly("degree", &bezier_t::degree_)
          .def_readonly("nbWaypoints", &bezier_t::size_)
      ;
  */
  /** END bezier curve**/

  /** BEGIN spline curve function**/
  class_<polynom_t>("polynomial",
                    init<const polynom_t::coeff_t, const real, const real>())
      .def("__init__", make_constructor(&wrapSplineConstructor))
      .def("min", &polynom_t::tmin)
      .def("max", &polynom_t::tmax)
      .def("__call__", &polynom_t::operator())
      .def("derivate", &polynom_t::derivate);
  /** END cubic function**/

  /** BEGIN spline curve**/
  class_<spline_t>("spline", no_init)
      .def("__init__", make_constructor(&wrapExactCubicConstructorvoid))
      .def("__init__", make_constructor(&wrapExactCubicConstructorPolySequence))
      .def("__init__", make_constructor(&spline_by_concatenation_constructor))
      .def("min", &spline_t::tmin,
           bp::return_value_policy<bp::return_by_value>())
      .def("max", &spline_t::tmax,
           bp::return_value_policy<bp::return_by_value>())
      .def("__call__", &spline_t::operator(),
           bp::return_value_policy<bp::return_by_value>())
      .def("derivate", &spline_t::derivate,
           bp::return_value_policy<bp::return_by_value>())
      .def("create_spline_from_waypoints", &spline_from_waypoints,
           boost::python::args("waypoints", "time vector"),
           "Creates a cubic spline from a set of way points")
      .def("create_spline_from_waypoints_constr", &spline_from_waypoints_constr,
           boost::python::args("waypoints", "time vector", "constraints"),
           "Creates a spline from a set of way points and constraints")
      .def("load_from_file", &spline_t::loadFromFile,
           boost::python::args("filename"), "Loads *this")
      .def("save_to_file", &spline_t::saveToFile,
           boost::python::args("filename"), "Saves *this");
  /** BEGIN force curve**/

  class_<force_t>("forcecurve", no_init)
      .def("__init__", make_constructor(&wrapForceCurveConstructorvoid))
      .def("__init__", make_constructor(&wrapForceCurveConstructorSplines))
      .def("min", &force_t::tmin,
           bp::return_value_policy<bp::return_by_value>())
      .def("max", &force_t::tmax,
           bp::return_value_policy<bp::return_by_value>())
      .def("__call__", &force_t::operator(),
           bp::return_value_policy<bp::return_by_value>())
      .def("derivate", &force_t::derivate,
           bp::return_value_policy<bp::return_by_value>())
      .def("set_motion_vector", &force_t::setMotionVector,
           boost::python::args("motionvector"),
           "sets motion vector for derivate")
      .def("load_from_file", &force_t::loadFromFile,
           boost::python::args("filename"), "Loads *this")
      .def("save_to_file", &force_t::saveToFile,
           boost::python::args("filename"), "Saves *this");

  /** END force curve**/

  /** END bezier curve**/

  /** BEGIN curve constraints**/
  class_<curve_constraints_t>("curve_constraints", init<>())
      .add_property("init_vel", &get_init_vel, &set_init_vel)
      .add_property("init_acc", &get_init_acc, &set_init_acc)
      .add_property("end_vel", &get_end_vel, &set_end_vel)
      .add_property("end_acc", &get_end_acc, &set_end_acc);
  /** END curve constraints**/
}
}  // namespace parametriccurves
